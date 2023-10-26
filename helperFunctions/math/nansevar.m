function [dataSevar] = nansevar(input,dim)
% This function calculates the standard error of the variance using a
% bootstrapping method.
%
% Inputs:
% - input: [? x ... nreps .. x ?] - some tensor with nreps as the size of
% the dimension you want to take the standard error of the variance over
% - dim: scalar, the dimension with nreps that you want to get the SEVar
% over
%
% Outputs:
% - dataSevar = [? x ... 1 .. x ?] - the standard error of the variance
% calculated for the given dimension
%
% Adam Smoulder, 8/10/22

nboots = 1000; % we'll do 10000 bootstraps and do the standard dev of their variance

% Permute to put our dim in question first
inputSize = size(input);
ndims = length(inputSize);
nreps = size(input,dim);
% allDims = 1:10;
% permDimOrder = [dim setdiff(1:ndims,dim)];
permDimOrder = 1:ndims;
permDimOrder(1) = dim;
permDimOrder(dim) = 1;
permInput = permute(input,permDimOrder);
permSize = size(permInput);

bootVals = nan([nboots permSize(2:end)]);
for b = 1:nboots
    curInds = randsample(nreps,nreps,true); % sample with replacement
    bootData = permInput(curInds,:,:,:,:,:,:); % if we have more than 6 dims, we're fucked
    bootVariance = nanvar(bootData);
    bootVals(b,:,:,:,:,:,:) = squeeze(bootVariance); 
end; clear b
depermOrder = 1:ndims;
depermOrder(1) = dim;
depermOrder(dim) = 1;
permSevar = std(bootVals);
dataSevar = permute(permSevar,depermOrder);
end

