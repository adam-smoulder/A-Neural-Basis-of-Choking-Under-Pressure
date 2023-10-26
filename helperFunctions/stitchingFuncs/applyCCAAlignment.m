function [trialData] = applyCCAAlignment(trialData,ccaStitchModel)
% This function takes a ccaStitchModel object (obtained through the script 
% STITCHING_usingCCAMethodAfterNormalStitching.m) and performs the
% alignment on the spike rates in trialData. It replaces "factorScores"
% with the CCA aligned projections.
%
% Inputs:
% - trialData: output structure of neural stitching
% - ccaStitchModel: structure output CCA alignment between days. Should
% have both the spike rate means and the alignment matrix for each day.
%
% Outputs:
% - trialData: same as input, but the "factorScores" field for each day has
% been replaced with the CCA aligned projections
%
% Adam Smoulder, 8/31/21

mu = ccaStitchModel.mu;
C = ccaStitchModel.c;
for i = 1:length(trialData)
    curSpikeMat = trialData(i).neuralData.spikeMatrix*1000; % convert to rate by using
    curDay = trialData(i).day;
    trialData(i).neuralData.factorScores = ((curSpikeMat-mu{curDay})'*C{curDay})';
end; clear i

end

