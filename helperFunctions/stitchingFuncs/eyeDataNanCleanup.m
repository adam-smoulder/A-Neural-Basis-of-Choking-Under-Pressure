function [eyeMatrix] = eyeDataNanCleanup(eyeMatrix)
% This function does basic outlier detection for eyeData and replaces nan
% values with a linear interpolation of the two surrounding points, then
% does a relatively coarse LPF on it. 
%
% Inputs:
% - eyeMatrix: [N x 2] matrix with X and Y eye data in each column
%
% Outputs:
% - eyeMatrix: [N x 2] matrix that is the cleaned up input
%
% Assumes eye matrix has already been centered and is oriented correctly.
%
% Adam Smoulder, 6/16/20

% Set any indices looking wayyyy off to nan in the first column
eyeMatrix(abs(eyeMatrix(:,1)) >= 200,1) = nan;
eyeMatrix(abs(eyeMatrix(:,2)) >= 200,1) = nan;

% Go through each index that has a nan in either X or Y and interpolate
nanInds = find(isnan(eyeMatrix(:,1)) | isnan(eyeMatrix(:,2)));
for k = 1:length(nanInds)
    firstInd = find(~isnan(eyeMatrix(1:nanInds(k)-1,1)),1,'last');
    lastInd = nanInds(k) + find(~isnan(eyeMatrix(nanInds(k)+1:end,1)),1,'first');
    if nanInds(k) == 1 % the very first index is bad...
       eyeMatrix(nanInds(k),:) = eyeMatrix(lastInd,:);
    elseif isempty(lastInd) % there's no non-nans after this
        eyeMatrix(nanInds(k),:) = eyeMatrix(firstInd,:);
    else % normal case
        eyeMatrix(nanInds(k),:) = 0.5*eyeMatrix(firstInd,:)+0.5*eyeMatrix(lastInd,:);
    end
end; clear k

% Smooth with 4th order forward/backward (8th order total) LPF @300Hz
Wn = 300/1000;                   
[b,a] = butter(4,Wn,'low'); 
eyeMatrix = filtfilt(b,a,eyeMatrix);
end

