function [goodnessScore,distFromBoundary,niterToComplete] = calculateGoodnessScore_trainTest(trainData,trainLabels,testData,testLabels)
% This function calculates "reach goodness" based on (1) an input data
% matrix of behavioral metrics, 1 per trial, and (2) trial status input
% (success, undershoot, overshoot). The general method is:
% - Make a QDA classifier to predict success label
% - For each trial, calculate distance from the success/oshoot and the
% success/ushoot decision boundary
%   - I couldn't find a closed form solution for this?? So I just use
%   gradient descent.
% - Sign the distance based on point location (+ if in success classifier
% zone, - if not)
% - Goodness score for a trial = minimum of the 2 distances
%
% Hypothetically this should work for any number of failure types, just so
% long as your first label is the one you want to compare to. 
%
% Inputs:
% - trainData: [ndata x ndims] matrix with the input metrics to train the
% classifier on. Recommended that this is balanced across extrinsic
% conditions that you aren't labeling (i.e., reward).
% - trainLabels: [ndata x 1] vector with input for each trial as 1 (success),
% 2 (undershoot), or 3 (overshoot) for the training data
% - testData: [ndata x ndims] data that you wish to test on to get out
% goodnessScores
% - testLabels: [ndata x 1] vector, same style as trainLabels but for
% testData.
%
% Outputs:
% - goodnessScore: [ndata x 1] score for test data
% - distFromBoundary: [ndata x 2] value with the (unsigned) distance from
% each of the decision boundaries for each point (first col is
% succ/undershoot boundary, second is succ/overshoot boundary)
% - niterToComplete: [ndata x 2] number of iterations for the given data
% point and boundary to calculate the distance
%
% Adam Smoulder, 7/16/21


[ndata,ndims] = size(testData);
classes = unique(trainLabels);
nclasses = length(classes);
discMdl = fitcdiscr(trainData,trainLabels,'discrimtype','quadratic','prior',ones(nclasses,1)/nclasses);

% Gradient descent parameters
nitermax = 1E7; % Max number of iterations for a given trial
thresh = 1E-6;  % Threshold for norm difference of output
allLR = [5E-3 4E-3 3E-3 1E-3 5E-4 1E-4 5E-5 1E-5 5E-6 1E-6]; % Learning rate; picked empirically
LR_TRANS = [10 2500 50000 500000 1E6 5E6];  % Number of iters after which to update learning rate

% Find distance to nearest point on the boundary for each class
distFromBoundary = nan(ndata,nclasses-1);
niterToComplete = nan(ndata,nclasses-1);
for c = 1:(nclasses-1)
    curCoeffs = discMdl.Coeffs(1,c+1);
    b0 = curCoeffs.Const; % scalar
    b1 = curCoeffs.Linear';  % 1 x ndims vector
    b2 = curCoeffs.Quadratic'; % ndims x ndims matrix (should be symmetric, transposing just for consistency)
    
    % Through some pen and paper work, we found an objective using
    % a Lagrange multiplier to find the point on the high-D parabola that
    % is closest to point xo:
    % L(x,lambda) = ||x-xo||^2-lambda*(b0+b1*x+x'*b2*x)
    % dL(X,lambda) = [2*(x-xo)+lambda*(b1+2*b2*x) ;  b0+b1*x+x'*b2*x ] = gradient, ndims+1 x 1
    % HL(X,lambda) = [[2*(I-lambda*b2)] b1+2*b2*x];     = hessian = ndims+1 x ndims+1
    %                [(b1+2*b2*x)'         0     ]]
    % So we'll do gradient descent on this
    lrMod = 0; % how many LR indices to skip; useful if exploding gradient bc too high LR
    i = 1; % because we have to sometimes redo an index if the LR is too high on a fit, we have to use a while loop
    while i <= ndata
        x0 = testData(i,:)';
        
        % Pick a guess for x and lambda; we assume x is a column vector
        lambda = 1;
        x = x0;
        
        % Iterate
        count = 1;
        x_prev = zeros(size(x));
        lambda_prev = 0;
        allXVals = nan(nitermax,ndims);
        allLambdaVals = nan(nitermax,1);
        allNormGradVals = nan(nitermax,1);
        while (count < nitermax) && (norm([x ; lambda]-[x_prev ; lambda_prev])>thresh)
            allXVals(count,:) = x;
            allLambdaVals(count) = lambda;
            lr = allLR(length(LR_TRANS)+1-sum(count<LR_TRANS)+lrMod); % Picks the correct learning rate based on trial count
            
            % Calculate gradient at the current guesses
%             L = sqrt(sum((x-x0).^2))-lambda*(b0+b1*x+x'*b2*x);
            dL = [2*(x-x0)+lambda*(b1'+2*b2*x) ; b0+b1*x+x'*b2*x];
            allNormGradVals(count) = norm(dL);
            
            % update count and vals
            x_prev = x;
            lambda_prev = lambda;
            x = x - lr*dL(1:ndims);
            lambda = lambda + lr*dL(end);
            count = count+1;
        end
        if count == nitermax
            disp(['Hit nitermax for i c = ' num2str(i) ' ' num2str(c) ', skipping (leaving nan)'])
            i = i+1;
            continue
        end
        
        % The distance from the boundary is L (since the lambda term should
        % effectively be 0). 
        distFromBoundary(i,c) = sqrt(sum((x-x0).^2))-lambda*(b0+b1*x+x'*b2*x);
        niterToComplete(i,c) = count;
        
        if isnan(distFromBoundary(i,c)) && lrMod < 3
            disp(['Nan hit, i c = ' num2str(i) ' ' num2str(c) ', trying lower LR'])
            i = i-1;
            lrMod = lrMod+1; % lower the learning rates
        elseif lrMod >=3
            disp(['Failed 3 attempts for i c = '  num2str(i) ' ' num2str(c) ', skipping (leaving nan)'])
            lrMod = 0;  % reset it
        else
            lrMod = 0; % reset it
        end
        
        % Update iterating counter
        if mod(i,100)==0, disp([c i]); end
        i = i+1;
    end; clear i
end; clear c

% Give the distance a sign based on whether it's in the success zone or not
meanSuccessPt = mean(testData(testLabels==1,:));
signedDistFromBoundary = nan(size(distFromBoundary));
for c = 1:nclasses-1
    curCoeffs = discMdl.Coeffs(1,c+1);
    b0 = curCoeffs.Const; % scalar
    b1 = curCoeffs.Linear';  % 1 x ndims vector
    b2 = curCoeffs.Quadratic'; % ndims x ndims matrix (should be symmetric, transposing just for consistency)
    x = meanSuccessPt';
    succSign = sign(b0+b1*x+x'*b2*x);
    
    x = testData';
    signedDistFromBoundary(:,c) = distFromBoundary(:,c).*sign(diag(b0+b1*x+x'*b2*x))*succSign;
end; clear c


% Finally, goodness score is the minimum of these. We take the exponential
% to make the distribution for successes more normal.
goodnessScore = exp(min(signedDistFromBoundary,[],2));
disp('Goodness score calculated!')

end

