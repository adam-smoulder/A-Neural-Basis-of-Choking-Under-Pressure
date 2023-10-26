function [FRzeroFactorMeans] = getFRZeroScoresFromStitchedModel(stitchingModel,dayLabels)
% This function trys to get at what a "zero" firing rate would look like in
% our factors after stitching. Because we de-mean individual unit firing
% rates each day, this isn't 100% trivial, but we do the following to try
% to keep it as fair/real as possible:
% 1) make a matrix of fake FRs, one per neuron/trial, that are zeros if a
% neuron is present on the day that trial was recorded and nan otherwise
% 2) subtract the actual means for each unit by each day (just like we do 
% for real data that we put through the stitching model)
% 3) put the data through the FA model -> orthogonalize -> take the mean
%
% What this basically means is that we have a weighted average of what
% zero-d out FR would be across days, where each day's weight is
% proportional to how many trials it has.
%
% Inputs:
% - stitchingModel: structure acquired after normal stitching preprocessing
% was done that has the stitching parameters from the actual toolbox, along
% with day means (d_day)
% - dayLabels: [ntrials x 1] array with a label for which day each trial
% belongs to
%
% Outputs:
% - FRzeroFactorMeans: [1 x nfactors] vector with the factor scores
% associated with FR = 0.
%
% Adam Smoulder, 1/27/21

nunits = length(stitchingModel.d);
ntrials = length(dayLabels);
days = unique(dayLabels);
ndays = length(days);
FRzeroData = zeros(ntrials,nunits);
dayMeans = stitchingModel.d_day;  % NOTE we already have nans in these means for when a unit isn't present!
for a = 1:ndays
    curInds = dayLabels==days(a);
    FRzeroData(curInds,:) = FRzeroData(curInds,:)-dayMeans(:,a)'; % subtracting by nans = nan, so units not present on a given day are nan'd
end; clear a

% Get FA scores
FRzeroFactorScores_raw = inferFALatents(stitchingModel, FRzeroData)';

% Put latents into an orthonormalized coordinate system and average
[~, FRzeroFactorScores] = orthonormalizeFAMdl(stitchingModel, FRzeroFactorScores_raw');
FRzeroFactorMeans = mean(FRzeroFactorScores);

end

