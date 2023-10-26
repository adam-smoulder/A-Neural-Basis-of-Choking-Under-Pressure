function [predLabels,LL_byClass] = gaussMLE(x,mu,sigma,varargin)
% For datapoint(s) x, evaluates which class of input has the highest
% log-likelihood.
%
% Inputs:
% - x: [npts x ndims] data to evaluate
% - mu: [nclasses x ndims] means of gaussian distributions
% - sigma: [nclasses x ndims x ndims] covariances of gaussian distribution
%
% Optional Inputs:
% - invSigma: [nclasses x ndims x ndims] inverse of covariance matrices. 
% Having this provided can speed things up
%
% Outputs:
% - predLabels: [npts x 1] selected class
% - LL_byClass: [npts x nclasses] log-likelihood for each class
%
% Adam Smoulder, 9/20/21

if nargin > 3
    invSigma = varargin{1};
else
    invSigma = nan(size(sigma));
    for c = 1:size(sigma,1)
        invSigma(c,:,:) = inv(squeeze(sigma(c,:,:)));
    end; clear c
end

N = size(x,1); % npts, ndims
C = size(mu,1); % nclasses

% Get multivariate gaussian LL for each class
LL_byClass = nan(N,C);
for c = 1:C
    LL_byClass(:,c) = gaussLL(x,mu(c,:),squeeze(sigma(c,:,:)),squeeze(invSigma(c,:,:)));
end; clear c

% Pick highest value
[~,predLabels] = max(LL_byClass,[],2);
predLabels(isnan(sum(LL_byClass,2))) = nan; % any nan LLs = make output nan.
end

