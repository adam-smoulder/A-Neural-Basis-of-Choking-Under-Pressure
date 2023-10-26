function [lowBound,highBound] = bootConfInt(data,percent)
% bootConfInt returns the confidence interval bounds corresponding to 
% whatever is in "percent" (e.g. for a 95% CI, percent = 95). This is
% assuming the input data is coming from a bootstrapped distribution
% (distribution of means).
%
% data is assumed to be a 1D vector of length nboots. data can also be a D 
% dimensional matrix, where this operation is then  performed along the 
% LAST dimension. Note that because this is assumed to come from 
% bootstrapped data, the input data is assumed to represent means, meaning 
% that the CI is read directly from them (not normalized by the length of 
% the data or anything)
%
% Inputs:
% - data: input data, assumed to be means from a bootstrapped distribution
% - percent: what percent confidence interval is desired. If data is a
% matrix, place the dimension with the bootstraps AT THE END (so if I was
% using data that was tubes x rewards x bootstraps, it should be of size
% either ntubes x nrewards x nboots or nrewards x ntubes x nboots)
%
% Outputs:
% - lowBound: lower bound of CI
% - highBound: upper bound of CI
% I suppose I could've just output a [2x1] vector or something...oh well
%
% Adam Smoulder, 3/15/19 (edit 4/17/19)

dims = size(data);
ndims = sum(dims~=1);

if ndims == 1 % vector input; output is scalar
    lowBound = squeeze(prctile(data,(100-percent)/2));
    highBound = squeeze(prctile(data,100-(100-percent)/2));
else % matrix input; output is matrix with the top D-1 dimensions of the input
    lowBound = squeeze(prctile(data,(100-percent)/2,ndims));
    highBound = squeeze(prctile(data,100-(100-percent)/2,ndims));
end 

end

