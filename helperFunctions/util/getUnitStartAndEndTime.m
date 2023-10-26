function [unitStartTime,unitEndTime] = getUnitStartAndEndTime(data)
% This function finds the first and last time bin each unit is present
% (meaning it is NOT nan) for the given input data matrix, which is assumed
% to be # of time points X # of units (same as stitchingInput.obs
% structure)
%
% Adam Smoulder, 6/6/19

[ntime, nunits] = size(data);
unitStartTime = nan(nunits,1);
unitEndTime = nan(nunits,1);
for i = 1:nunits
    unitStartTime(i) = find(~isnan(data(:,i)),1,'first');
    firstNanTime = find(isnan(data(unitStartTime(i):end,i)),1,'first')+unitStartTime(i)-1;
    if isempty(firstNanTime) % it lasts until the end
        unitEndTime(i) = ntime;
    else % it stops at some point
        unitEndTime(i) = firstNanTime-1;
    end
end; clear i
end

