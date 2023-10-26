function [theMSE] = mseBAD(predData,realData)
% Just calculates mean squared error. Nothing fancy. MATLAB's one toolbox
% seems to have a built in version with the same name, but it seems to have
% a buncha bells and whistles that I really don't give a shit about.
%
% This function only calculates MSE over the first dimension of whatever
% inputs you use (which are assumed to be the same size). If either the
% prediction or the real data are NaN for a given point, it's ignored for
% the calculation.
%
% Inputs:
% - predData: N x ??? vector/matrix with the prediction of the data across
% N datapts
% - realData: N x ??? vector/matrix with the real data values
%
% Outputs:
% - theMSE: scalar (or matrix of the same size as the input after the first
% dimension) with the mean-squared error
%
%
% Adam Smoulder, 7/10/19

% theMSE = nanmean((predData-realData).^2); % woooo one whole line!
% theMSE = nanmean(abs(predData-realData)); % if for some odd reason you wanted a command called mse to get mean absolute error...

% or if you want r-square...
SSE = sum((realData-predData).^2);
SST = sum((realData-mean(realData)).^2);
theMSE = 1-SSE./SST;
disp("USING R^2")
end

