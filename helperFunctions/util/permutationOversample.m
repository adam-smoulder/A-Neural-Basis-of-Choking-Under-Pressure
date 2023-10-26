function [OSInds] = permutationOversample(n,n_OS)
% This function oversamples values on 1:n up to n_OS values. As opposed to
% voersampling randomly, we oversample by just stacking random permutations
% of n until we get more datapts than n_OS, then we chop off the excess.
%
% The main utility of this is for data balancing for parametric models.
%
% Inputs:
% - n: scalar, number of original samples
% - n_OS: scalar, number of values to oversample to
%
% Output:
% - OSInds: [n_OS x 1] oversampled indices
%
% For instance, if you pass in permutationOversample(3,8), you may get
% [1 2 3 2 3 1 1 3]'. 1 2 3 is repeated (randomized) twice, then a truncated
% version gets up to n_OS
%
% Adam Smoulder, 7/21/22

nreps = ceil(n_OS/n);
OSInds = nan(nreps*n,1);
for rep = 1:nreps
    OSInds((1:n)+(rep-1)*n) = randperm(n);
end; clear rep
OSInds = OSInds(1:n_OS);
end

