function [validUnitInds,nminByUnit] = unitsThatHaveMinTrialCount(data,nmin,labels)
% This function identifies which units have at least nmin trials for each
% label condition.
%
% Inputs:
% - data: [ntrials x nunits] data matrix
% - nmin: scalar, number of minimum trials for each label condition
% - labels: [ntrials x nlabelcols] labels for each trial. If nlabels > 1, each
% unique combination of labels will be considered (i.e., if you have
% [directionLabels rewardLabels], we would look at dir==1&rew==1,
% dir==1&rew==2, etc.)
%
% Outputs:
% - validUnitInds: [1 x nunits] boolean array with the valid units
% identified
% - nminByUnit: [1 x nunits] array with the number of trials in the minimum
% condition present for this unit
%
% Adam Smoulder, 8/2/22

% Combine labels; assumes positive integers
if size(labels,2) > 1
    nlabelcols = size(labels,2);
    maxVals = max(labels);
    labels = labels*[1 cumsum(maxVals(1:end-1))]'-maxVals(1);
end

% Make labels ordered 1 to max unique val for ease, then count
[~,labels] = ismember(labels,unique(labels));
nlabels = length(unique(labels));

% For each unit, identify how many trials of each label are present, and
% store the minimum
[~,nunits] = size(data);
nminByUnit = nan(nunits,1);
for u = 1:nunits
    curData = data(:,u);
    validInds = ~isnan(curData);
    curLabelCounts = histcounts(labels(validInds),0.5:nlabels+0.5);
    nminByUnit(u) = min(curLabelCounts);
end; clear u

% Identify if the min is above the input value
validUnitInds = nminByUnit >= nmin;

end

