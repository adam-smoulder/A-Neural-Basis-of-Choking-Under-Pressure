function [allBoots, bootMean, boot95low, boot95up, bar95low, bar95up, SE] = bootstrapBinaryEvent(eventCounts,totalCounts,nboots)
% This function is for bootstrapping a binary (bernoulli) event, where you
% have counts of how often this event happened along with the total number
% of trials.
%
% Inputs:
% - eventCounts: matrix with the counts of how often the event happened
% - totalCounts: matrix (same size as eventCounts) with counts of how many
% opportunities there were for the event to happen
% - nboots: how many bootstraps to make
%
% Outputs:
% - allBoots: [nboots x size(eventCounts)] matrix with all bootstraps'
% counts
% - bootMean: [size(eventCounts)] matrix with the average of the bootstraps
% - boot95low: 2.5 percentile of results for 95% CI (same size as bootMean)
% - boot95up: 97.5 percentile of results for 95% CI (same size as bootMean)
% - bar95low: how far below the mean the 2.5 percentile is
% - bar95up: how far above the mean the 97.5 percentile is
% = SE: standard error of the metric (stdev of the bootstrap)
%
% These last 2 are redundant with the 2 before it; I just use em bc that's
% what MATLAB's errorbar function takes.
%
% Adam Smoulder, 6/2/20 (edit 7/27/20)

% We're going to just go element by element for doing this and then reshape
% stuff properly at the end. 
dataShape = size(eventCounts);
eventCounts = eventCounts(:);
totalCounts = totalCounts(:);
nelements = length(eventCounts);
bootCounts = nan(nboots,nelements);

% For each element, do the bootstrap
for e = 1:nelements
    for b = 1:nboots
        % sample WITH replacement; if the "index" is above the # of events, we say it didn't happen. 
        % This is the same as labeling totalCounts(e) with either a 1 or 0, where the sum is eventCounts(e), then redrawing
        % ...which is the same as doing a Bernoulli draw with p = eventCounts(e)/totalCounts(e) a total of totalCounts(e) times
        bootCounts(b,e) = sum(randsample(totalCounts(e),totalCounts(e),true)<=eventCounts(e)); 
    end; clear b                                                             
end; clear e                                                                 

% Reshape into the proper form for output and make outputs
allBoots = reshape(bootCounts,[nboots dataShape]);
bootMean = squeeze(mean(allBoots));
boot95low = squeeze(prctile(allBoots,2.5));
boot95up = squeeze(prctile(allBoots,97.5));
bar95low = bootMean-boot95low;
bar95up = boot95up-bootMean;
SE = squeeze(std(allBoots));
end

