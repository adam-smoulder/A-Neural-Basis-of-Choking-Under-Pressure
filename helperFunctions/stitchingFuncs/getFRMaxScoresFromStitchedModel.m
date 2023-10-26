function [FRmaxUnitMeans] = getFRMaxScoresFromStitchedModel(taskInfo,trialData)
% This function gets the max firing rate for each unit for each day and
% applies the stitchingModel to it.
%
% Adam Smoulder, 1/27/21

stitchingModel = taskInfo.stitchingModel;
dayLabels = [trialData.day]';
days = unique(dayLabels);
ndays = length(days);
ntrials = length(trialData);
nunits = length(stitchingModel.d);
dt = median(diff(trialData(1).time)); % in ms
binLength = round(200/dt); % 200 ms
kernel = ones(1,binLength)/binLength*1000; % convert to Hz

% Get the max FR for each unit for each day
FRmaxUnitMeans = nan(ntrials,nunits);
for a = 1:ndays
    curInds = dayLabels==days(a);
    maxFRMat = [];
    curData = trialData(curInds);
    for i = 1:length(curData)
        trialFRMat = convByDim(curData(i).neuralData.spikeMatrix,kernel',2,'same');
        maxFRMat = max([maxFRMat trialFRMat],[],2);
    end; clear i
    FRmaxUnitMeans(curInds,taskInfo.unitIDs{a}) = repmat(maxFRMat',[sum(curInds),1]);
end; clear a

% Subtract day-means
dayMeans = stitchingModel.d_day;  % NOTE we already have nans in these means for when a unit isn't present!
for a = 1:ndays
    curInds = dayLabels==days(a);
    FRmaxUnitMeans(curInds,:) = FRmaxUnitMeans(curInds,:)-dayMeans(:,a)';
end; clear a

% Get FA scores
FRmaxFactorScores_raw = inferFALatents(stitchingModel, FRmaxUnitMeans)';

% Put latents into an orthonormalized coordinate system and average
[~, FRmaxFactorScores] = orthonormalizeFAMdl(stitchingModel, FRmaxFactorScores_raw');
FRmaxUnitMeans = mean(FRmaxFactorScores);
end

