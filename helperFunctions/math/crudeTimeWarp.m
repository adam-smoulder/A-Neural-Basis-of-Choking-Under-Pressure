function [newTrace] = crudeTimeWarp(baseTrace,twf,sf)
% This function does a REALLLLY shitty linear time warp and scaling. REALLY
% shitty. No averaging, lots of rounding, it's called "crudeTimeWarp" not
% "wellImplementedAndWellThoughtThroughTimeWarp".
%
% Inputs:
% - baseTrace: [nbins x 1] trace that you want to warp
% - twf: scalar time-warp factor. >1 = stretching, <1 = compressing.
% - sf: scalar scaling factor
%
% Outputs:
% - newTrace: [nbins x 1] vector with your time-warped and scaled trace.
%
% Adam Smoulder, 1/18/21

nbins = length(baseTrace);
oldInds = 1:nbins;
newInds = round(1:twf:nbins);
newTrace = nan(size(baseTrace));

if twf >= 1 % stretching & linear interpolation between new pts.
    oldInds(length(newInds)+1:end) = [];
    newTrace(newInds) = sf*baseTrace(oldInds);
    nanInds = find(isnan(newTrace));
    newTrace(nanInds) = interp1(newInds,newTrace(newInds),nanInds);
else % compressing & repeating last value for the end
    newInds(length(oldInds)+1:end) = [];
    [~,indsWithinNewInds] = unique(newInds);
    newTrace(newInds(indsWithinNewInds)) = sf*baseTrace(oldInds(indsWithinNewInds));
    newTrace(max(newInds)+1:end) = baseTrace(oldInds(end));
end

end

