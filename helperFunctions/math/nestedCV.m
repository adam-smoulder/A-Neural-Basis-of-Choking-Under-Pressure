function [estMean,estSEM,normalCVMean,minVal,maxVal] = nestedCV(X,Y,method,K,R,varargin)
% This function performs nested cross-validation per Bates, Hastie, and
% Tibshirani 2021. The main goal is to get a confidence interval (and
% hence, standard error bars) for the result of some cross-validated model
% fitting.
%
% The general gist is that the common procedure to calculate error bars as
% S.E. over folds underestimates the actual variance of the statistic
% because it treats the fold-results as independent; they are not
% independent, as the data used to fit the test data in each fold overlaps,
% and the test data from one fold are used in the training model for
% another. To circumvent this, the paper shows that the mean squared error
% of a cross-validated metric based on specific training data can be
% estimated by performing cross-validation K-1 folds and comparing results
% with the held out fold (over many repetitions).
%
% Honestly, the paper was kind of dense for me and took me a good 3 long
% readthroughs to understand. I'll try to explain the algorithm steps in 
% the comments as we go. I'll also try to use the same variable conventions
% as the paper when possible.
%
% One difference: I balance the data by randomly subsampling here. While
% this does lose some of the immediate guarantees of the equations in the
% paper, this should yield a conservative estimate of error (greater than
% actual).
%
% Inputs:
% - X: [N x D] data, where N is number of observations and D is number of
% dimensions. 
% - Y: [N x 1] labels that we are decoding form the data.
% - method: string, the classification/regression method to be used.
% Options include 'GNB','LDA'
% - K: scalar, number of folds to use (must be >= 3)
% - R: scalar, number of repetitions to use.
%
% Optional inputs:
% - LOGGING: boolean; if true, will show update every 50 reps (default:
% true)
%
% Outputs:
% - estMean: scalar, the mean quantity of interest. For 'GNB' and 'LDA',
% this is decoding accuracy (%).
% - estSEM: scalar, the S.E. of the mean quantity of interest.
%   - Importantly, when comparing models of different #s of trials, the
%   estSEM from this is valid, but you should separately run
%   cross-validation on balanced data over many reps to avoid bias in some
%   models due to higher trial counts.
% - normalCVMean: scalar, output from normal cross-validation
% - minVal: minimum value for error bars; as if you had N independent
% samples
% - maxVal: max value for error bars; as if you had N/K independent samples
%
% Algorithm source: Bates, Hastie, and Tibshirani, "Cross-validation: what
% does it estimate and how well does it do it?", arXiv, 2021 April 15
%
% Adam Smoulder, 10/4/21


if ~isempty(varargin)
    LOGGING = varargin{1};
else
    LOGGING = true;
end


% Get overall size of data and group data
[N,D] = size(X);
[~,inds_byLabel] = groupDataByLabel(X,Y);
Yvals = unique(Y); % needed in case Y values are not 1:Cy
n_byLabel = cellfun(@(x) length(x), inds_byLabel);
n_mincond = min(n_byLabel(:)); % number of points to subsample to

% Get info about labels
classes = unique(Y);
C = length(classes); % number of classes

% We will scramble our indices each time, so we only have to generate fold
% labels once
foldLabels = repmat(crossvalind('kfold',n_mincond,K),[C 1]);

% Run through each repetition
es = nan(R,K,C*ceil((K-1)/K*n_mincond)); % hold out error; each point is tested once per rep
a_list = nan(R,K); % First term of MSE, difference between inner-CV estimate and hold out error (squared)
b_list = nan(R,K); % Second term of MSE, basically the squared standard error of the hold out error for a given rep x fold
e_out_all = nan(R,N); % all outer CV test errors ; this is the same as results from normal CVing, meaning one label per point
for r = 1:R % for each rep
    % Scramble our indices within each class and subsample this rep's data
    repPerm = cellfun(@(x) randperm(size(x,1))', squeeze(inds_byLabel),'UniformOutput',false);
    repInds_cell = cellfun(@(x,y) x(y(1:n_mincond),:), squeeze(inds_byLabel), repPerm, 'UniformOutput', false);
    repInds = cellArrayToVector(repInds_cell);
    
    % Now do nested cross-validation
    for k = 1:K % outer CV Loop
        innerFolds = setdiff(1:K,k);
        e_in = nan(n_mincond*C,1); % should be C-1, but indexing with that is way too hard, so we'll just remove nans at the end
        for ki = innerFolds % inner CV Loop
            % Get train/test inds together
            repTrainInds = ~ismember(foldLabels,[k ki]); % indices in our subsample
            trainInds = repInds(repTrainInds);  % indices in the full data
            repTestInds = foldLabels==ki;
            testInds = repInds(repTestInds); 
            testY = Y(testInds);
            
            % Train and test for this inner fold
            switch method
                case 'GNB'
                    [modelMeans,~,modelDiagCovs] = fitGNBModel(X(trainInds,:), Y(trainInds)); % train
                    predY = gaussMLE(X(testInds,:), modelMeans, modelDiagCovs);
                case 'LDA'
                    [~, ~, ~, ~,modelMeans,modelNoiseCov] = linearDiscAnalysis(X(trainInds,:), Y(trainInds),'DO_CV',false);
                    predY = gaussMLE(X(testInds,:), modelMeans, repmat(permute(modelNoiseCov,[3 1 2]),[C 1 1]));
                otherwise
                    error('Method not recognized')
            end
            e_in(repTestInds) = 1-(Yvals(predY)==testY); % 1/0 error
        end; clear ki
        e_in(isnan(e_in)) = []; % outer test fold has no values in this; remove them.
        
        
        % Get training and testing inds for the outer CV
        repTrainInds = ~ismember(foldLabels,k); % indices in our subsample
        trainInds = repInds(repTrainInds);  % indices in the full data
        repTestInds = foldLabels==k;
        testInds = repInds(repTestInds);
        testY = Y(testInds);
        
        % Fit a model on all inner-fold data and test on the outer fold
        switch method
            case 'GNB'
                [modelMeans,~,modelDiagCovs] = fitGNBModel(X(trainInds,:), Y(trainInds)); % train
                predY = gaussMLE(X(testInds,:), modelMeans, modelDiagCovs);
            case 'LDA'
                [~, ~, ~, ~,modelMeans,modelNoiseCov] = linearDiscAnalysis(X(trainInds,:), Y(trainInds),'DO_CV',false);
                predY = gaussMLE(X(testInds,:), modelMeans, repmat(permute(modelNoiseCov,[3 1 2]),[C 1 1]));
            otherwise
                error('Method not recognized')
        end
        e_out = 1-(Yvals(predY)==testY);
        
        % Get the values we need out of this 
        a_list(r,k) = (mean(e_in)-mean(e_out))^2;
        b_list(r,k) = var(e_out)/length(e_out); % Note the paper has an error and just says to use var(e_out); this is wrong (per the github code and eq 15)
        es(r,k,1:length(e_in)) = e_in; % note - if N/K is not an integer, then this size will differ across folds, leaving nans in this matrix 
        e_out_all(r,testInds) = e_out;
    end; clear k
    if mod(r,50)==0 && LOGGING, disp(['Completed rep ' num2str(r)]); end
end; clear r

% Get average outputs and get bias
MSE_est = mean(a_list(:))-mean(b_list(:));
err_NCV_est = nanmean(es(:)); % nanmean in the case of unequal fold sizes
err_CV_est = nanmean(e_out_all(:));
bias = (1+(K-2)/K)*(err_NCV_est-err_CV_est);

% And calculate what we need! We'll also change it from error to decode
% accuracy and convert it from fraction to percentage
estMean = err_NCV_est-bias;
estMean = 100*(1-estMean);
estSEM = sqrt((K-1)/K*MSE_est);
estSEM = 100*estSEM;
normalCVMean = (1-err_CV_est)*100;

% Finally, they add that the error should fall between a minimum and
% maximum range of the SE over the test data (a lower bound) and the SE
% over the test data as if you only have your # of samples in each fold
% (upper bound). So we apply that to our estSEM
minVal = 100*mean(nanstd(e_out_all,[],2))/sqrt(n_mincond*K);
maxVal = 100*mean(nanstd(e_out_all,[],2))/sqrt(n_mincond);
if estSEM < minVal
    disp('CORRECTING USING MINVAL (SE for N indep points)')
    estSEM = minVal;
elseif estSEM > maxVal
    disp('CORRECTING USING MAXVAL (SE for N/K indep points)')
    estSEM = minVal;
end

end

