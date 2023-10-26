function [stitchingInput,indsToKeep] = stitchingSelection_day(stitchingInput, daysToKeep)
% stitchingSelection methods are made to crunch down the data input to the
% stitching alignment algorithm by various labels. This one allows you to
% specify which days to keep.
%
% inputs:
% stitchingInput - structure with many fields pertinent to stitching and
% data labels
% daysToKeep - which days of data should be kept; to remove, simply instead
% input setdiff(1:ndays, daysToRemove) as the argument.
%
% Adam Smoulder, 1/10/19


if ~isempty(stitchingInput.masterInfo)
    indsToKeep = ismember(stitchingInput.masterInfo(:,1),daysToKeep);
    stitchingInput.masterInfo = stitchingInput.masterInfo(indsToKeep,:);
    stitchingInput.ntimebins = length(stitchingInput.masterInfo);
else
    indsToKeep = ismember(stitchingInput.assignments(:,6),daysToKeep);
    stitchingInput.assignments = stitchingInput.assignments(indsToKeep,:);
    stitchingInput.ntimebins = length(stitchingInput.assignments);
end
stitchingInput.lastTimeBins = stitchingInput.lastTimeBins(daysToKeep);
stitchingInput.filenames = stitchingInput.filenames(daysToKeep);
stitchingInput.lastTrials = stitchingInput.lastTrials(daysToKeep);
stitchingInput.spikeRates = stitchingInput.spikeRates(daysToKeep);
unitsToKeep = sum(~isnan(stitchingInput.obs))>0;
stitchingInput.obs = stitchingInput.obs(indsToKeep,unitsToKeep);
end

