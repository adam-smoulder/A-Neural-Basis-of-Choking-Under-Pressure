function [dataSemed] = nansemed(input,dim)
% This function calculates the standard error of the median using a
% bootstrapping method.
%
% Inputs:
% - input: [? x ... nreps .. x ?] - some tensor with nreps as the size of
% the dimension you want to take the standard error of the median over
% - dim: scalar, the dimension with nreps that you want to get the SEMed
% over
%
% Outputs:
% - dataSemed = [? x ... 1 .. x ?] - the standard error of the median
% calculated for the given dimension
%
% Adam Smoulder, 2/19/21

if isempty(input) || sum(~isnan(input(:)))==0
    dataSemed = nan;
    disp('Empty or all nans')
    return
end

nboots = 1000; % we'll do 10000 bootstraps and do the standard dev of their medians

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
    bootMedian = nanmedian(bootData);
    bootVals(b,:,:,:,:,:,:) = squeeze(bootMedian); 
end; clear b
depermOrder = 1:ndims;
depermOrder(1) = dim;
depermOrder(dim) = 1;
permSemed = std(bootVals);
dataSemed = permute(permSemed,depermOrder);
end

