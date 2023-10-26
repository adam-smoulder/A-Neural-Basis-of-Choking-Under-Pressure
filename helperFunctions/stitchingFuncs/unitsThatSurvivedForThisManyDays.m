function [units] = unitsThatSurvivedForThisManyDays(data,numDays)
% unitsThatSurvivedForThisManyDays takes in the stitchingInput.obs matrix
% (which is # of timebins x # of units)as "data" and a scalar "numDays".
% The output is the indices of the units that are present for at least the
% number of days indicated by numDays. 
% 
% Ex: I want all the units that were present for at least 5 days. I'd use:
% unitsThatSurvivedForThisManyDays(stitchingInput.obs,5)
%
% Adam Smoulder, 6/6/19

% Find unit start and end times
[ntime, nunits] = size(data);
[unitStartTimes, unitEndTimes] = getUnitStartAndEndTime(data);

% So long as there is at least one new unit per day, we can find the start
% and end time point for each day super quickly:
dayStartTimes = unique(unitStartTimes);
dayEndTimes = unique(unitEndTimes);

% Now for start and end days, we just find the indices for each unit of the
% day start and end times! Following this, # days present is easy
[~, unitStartDay] = ismember(unitStartTimes,dayStartTimes);
[~, unitEndDay] = ismember(unitEndTimes,dayEndTimes);
numDaysPresent = unitEndDay-unitStartDay+1;

% Finally, return only the units who have numDaysPresent >= numDays
units = find(numDaysPresent >= numDays);

end

