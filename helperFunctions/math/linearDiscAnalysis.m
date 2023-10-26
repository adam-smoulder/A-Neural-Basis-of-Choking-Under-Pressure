% A simple wrapper for MATLAB's linear discriminant analysis functionality
% that allows for lower numbers of dimensions to be specified.
%
% Inputs: 
%
%   input - a matrix of feature vectors.  Each row is a sample.  Should by 
%   S by D, where S is the number of samples and D is the dimensionality of
%   the data. 
%
%   labels - a vector of class labels.  Should be of length S. 
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format. 
%
%   NUMDIMOUT - Number of desired output dimensions (hacked in by Adam)
%   DO_CV - boolean on whether or not to do cross-validated decoding
%
% Outputs: 
%
%   w - A D by C-1 matrix, where C is the number of classes giving the
%   projection vectors. If the dimensionality of the data is less than C,
%   then C-1 will equal the dimensionality of the original data. 
%
%   objVal - value of the discriminant objective for the optimal projection
%
%   eigVls - Eigenvalues from the generalized eigenvalue problem on the
%   scatter matrices. Roughly represents the "importance" of each
%   projection dimension
%
%   cvDecodeAccs - simple 10 fold CV decoding accuracy. Doesn't try to
%   balance or anything
%
%   classMeans - A C x D matrix of the means for each class.
% 
%   sW - D x D "within class scatter" (or noise covariance). Average of the
%   within-condition covariances.
%
% NOTE with the way this is written, we assume uniform prior (and hence,
% don't use weighted average for calculating noise covariance)
%
% Adam Smoulder, 12/12/18 (edit 9/24/21)
%
function [w, objVal, eigVls, cvDecodeAccs, classMeans, sW] = linearDiscAnalysis(input, labels, varargin)

if size(input,1)~=size(labels,1)
    error(['Diff input and label sizes'])
end

NO_ASSERT = false; 
NUMDIMOUT = 0; % if 0, will do default
DO_CV = true; % do CVing for decode accuracy; defaults true

warnOpts(assignOpts(varargin)); 

% Get info on inputs
labelVals = unique(labels);
nlabels = length(labelVals);
[~,ndims] = size(input);

% % Method 1: Use MATLAB's fitting to get the covariance matrices
% discMdl = fitcdiscr(input,labels,'prior',ones(length(labelVals),1)/length(labelVals));
% sB = discMdl.BetweenSigma;
% sW = discMdl.Sigma;
% % % Importantly...at least per wikipedia, we should be normalizing our sB
% % by 1/C, not 1/(C-1) like MATLAB's implementation does. The model found
% % should be the same, but the eigenvalues and objective will be different.

% Method 2: Do it ourselves
classMeans = nan(nlabels,ndims);
classCovs = nan(nlabels,ndims,ndims);
for c = 1:nlabels
    curInds = labels == labelVals(c);
    classMeans(c,:) = mean(input(curInds,:));
    classCovs(c,:,:) = cov(input(curInds,:));
end; clear c
sW = squeeze(mean(classCovs));
sB = cov(classMeans,1);

% Get eigenvalues and fix order
[eigVecs, eigVls] = eig(sB, sW);
[eigVls, order] = sort(diag(eigVls),'descend');
eigVecs = eigVecs(:,order);

% find the rank (# labels - 1 or #dims) and extract w accordingly
m = min(length(unique(labels)), size(input,2));
eigVls = eigVls(1:(m-1));
w = eigVecs(:, 1:(m-1));
if NUMDIMOUT ~=0 % use the number of dimensions input
    w = w(:,1:NUMDIMOUT);
end

objVal = det(w'*sB*w)/det(w'*sW*w);


% Get 10 fold cross-validated decode accuracy
labelCounts = histcounts(labels,[labelVals-realmin ; max(labelVals)+realmin]);

% Make folds
nfolds = min(10,min(labelCounts));
testFoldInds = cell(nfolds,length(labelVals));
for i = 1:nlabels
    curLabel = labelVals(i);
    curInds = find(labels==curLabel);
    curInds = curInds(randperm(length(curInds))); % randomize order
    c = cvpartition(length(curInds),'KFold',nfolds);
    count = 0;
    for f = 1:nfolds
        testFoldInds{f,i} = curInds((count+1):(count+c.TestSize(f)));
        count = count+c.TestSize(f);
    end; clear f
end

% Test
cvDecodeAccs = nan(nfolds,1);
if DO_CV
    for f = 1:nfolds
        testInds = cellArrayToVector(testFoldInds(f,:));
        trainInds = setdiff((1:length(labels))',testInds);
        discMdl = fitcdiscr(input(trainInds,:),labels(trainInds),'prior',ones(length(labelVals),1)/length(labelVals));
%         discMdl = fitcdiscr(input(trainInds,:),labels(trainInds),'prior',ones(length(labelVals),1)/length(labelVals),'discrimtype','quadratic');
        testPredLabels = predict(discMdl, input(testInds,:));
        testRealLabels = labels(testInds);
        cvDecodeAccs(f) = 100*sum(testPredLabels==testRealLabels)./length(testPredLabels);
    end; clear f
end

end

