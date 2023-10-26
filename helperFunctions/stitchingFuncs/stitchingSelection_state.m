function [stitchingInput,indsToKeep] = stitchingSelection_state(stitchingInput, statesToKeep)
% stitchingSelection methods are made to crunch down the data input to the
% stitching alignment algorithm by various labels. This one allows you to
% specify which states to keep.
%
% inputs:
% stitchingInput - structure with many fields pertinent to stitching and
% data labels
% statesToKeep - which days of data should be kept; to remove, simply instead
% input setdiff(1:15, statesToRemove) as the argument.
%
% Adam Smoulder, 1/10/19

if ~isempty(stitchingInput.masterInfo)
    indsToKeep = ismember(stitchingInput.masterInfo(:,3),statesToKeep);
    stitchingInput.masterInfo = stitchingInput.masterInfo(indsToKeep,:);
    stitchingInput.ntimebins = length(stitchingInput.masterInfo);
    stitchingInput.obs = stitchingInput.obs(indsToKeep,:);
else
    % return the same thing - full trial doesn't have states
end
end

