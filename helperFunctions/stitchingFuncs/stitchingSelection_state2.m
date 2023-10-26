function [stitchingInput,indsToKeep] = stitchingSelection_state2(stitchingInput, statesToKeep_1targ, statesToKeep_2targ)
% stitchingSelection methods are made to crunch down the data input to the
% stitching alignment algorithm by various labels. This one allows you to
% specify which states to keep.
%
% inputs:
% stitchingInput - structure with many fields pertinent to stitching and
% data labels
% statesToKeep - which days of data should be kept; to remove, simply instead
% input setdiff(1:15, statesToRemove) as the argument.
% first arg is for 1 targ trials, 2nd ard is for 2 targ trials
%
% Adam Smoulder, 5/6/19

if ~isempty(stitchingInput.masterInfo)
    twoTargBins = stitchingInput.masterInfo(:,9);
    oneTargBins = ~twoTargBins;
    successBins = stitchingInput.masterInfo(:,8); % we only want to rout out bad states for successes; we keep failures
    statesToKeepBins_1targ = ismember(stitchingInput.masterInfo(:,3),statesToKeep_1targ);
    statesToKeepBins_2targ = ismember(stitchingInput.masterInfo(:,3),statesToKeep_2targ);
    indsToKeep = ((statesToKeepBins_1targ & oneTargBins) | (statesToKeepBins_2targ & twoTargBins)) | ~successBins;
    
    stitchingInput.masterInfo = stitchingInput.masterInfo(indsToKeep,:);
    stitchingInput.ntimebins = length(stitchingInput.masterInfo);
    stitchingInput.obs = stitchingInput.obs(indsToKeep,:);
else
    % return the same thing - full trial doesn't have states
end
end

