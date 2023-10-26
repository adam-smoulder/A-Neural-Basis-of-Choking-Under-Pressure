function [R2] = calculateR2(actualData,predictedData)
% This function calculates R^2 from actual data and predicted data. This
% only calculates R2 along the first dimension of the inputs; if you want
% to pool over other dimensions, you can always pass in actualData(:) or
% things like that.
%
% Inputs:
% - actualData: [nobs x ???] matrix with the actual data observations. Any
% number of dimensions can follow.
% - predictedData: [nobs x ???] matrix with predicted observations. Should
% be the same size as actualData.
%
% Outputs:
% - R2: [1 x ???] matrix with R2 for each dimension.
%
% Note: we calcualte R^2 the "correct" way where the maximum value is 1
% and negative values can be attained (if the variance of your residuals
% exceeds that of the original data - this means your prediction is
% probably VERY bad). Please don't do stuff like just squaring the
% correlation coefficient. That's bad news bears.
%
% Adam Smoulder, 1/18/21
SST = sum((actualData-mean(actualData)).^2); % total sum of squares
SSR = sum((actualData-predictedData).^2); % residual sum of squares
R2 = 1-SSR/SST;

end

