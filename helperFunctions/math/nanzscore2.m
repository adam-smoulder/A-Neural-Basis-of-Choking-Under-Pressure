function [output] = nanzscore2(data,varargin)
% This function does a median z-score of the data. We do this zscoring by
% using the mean instead of the median and using 0.741*iqr instead of the
% std. This is robust to outliers, though the median does have less sample
% power or whatever it's called (it takes more samples to get a good
% estimate than mean). Also ignores nans.
%
% Adam Smoulder, 7/12/22

if nargin==2
    dim = varargin{1};
else
    dim = 1;
end
output = (data-nanmedian(data,dim))./(0.741*iqr(data,dim));
end

