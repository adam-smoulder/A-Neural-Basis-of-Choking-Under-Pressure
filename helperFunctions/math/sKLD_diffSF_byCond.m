function [sKLD_byCond,sKLD_byNBins_byCond] = sKLD_diffSF_byCond(input,condLabels,successLabels,nbinsToTry)
% This function estimates the symmetric KL Divergence between the success
% and failure labels within each condition based on the given input metric.
%
% Inputs:
% - input: [ntrials x 1] input metric (i.e., reaction time, firing rate)
% - condLabels: [ntrials x 1] label that you want the data to be split by
% in evaluating sKLD. We'll balance counts for this.
% - successLabels: [ntrials x 1] label of success versus failure (what
% we'll be evaluating sKLD between)
% - nbinsToTry: a vector with how many bins to evaluate sKLD over. We
% average over these for the final output
%
% Outputs:
% - sKLD_byCond: [1 x nconds] vector of the sKLD between S/F of the given
% metric. We take the median over all nbins tried
% - sKLD_byNBins_byCond: Same thing but shown for all nbins tried
%
%
% Adam Smoulder, 6/15/22

% % Get inds for oversampling to balance conditions (Old method; produces
% % inconsistent results, though over many reps I guess would be
% % statistically better)
% [~,inds_byCond,n_byCond] = groupDataByLabel(input,condLabels);
% nconds = length(n_byCond);
% nOS = 10*max(n_byCond);
% OSInds = cellArrayToVector(cellfun(@(x) x(randi(length(x),[nOS 1])), inds_byCond,'uniformoutput',false));

% To balance conditions, we multiple count data from less common conditions (as best as we can)
conds = unique(condLabels); nconds = length(conds);
[~,~,n_byCond] = groupDataByLabel(input,condLabels);
nrecount_byCond = floor(max(n_byCond)./n_byCond);
balancedInput = input;
for c = 1:nconds
    if nrecount_byCond(c)~=1
        balancedInput = [balancedInput ; repmat(input(condLabels==conds(c)),[nrecount_byCond(c)-1 1])]; 
    end
end; clear c

% Get sKLD for each # of bins to try, then average over that
nNBinsToTry = length(nbinsToTry);
sKLD_byNBins_byCond = nan(nNBinsToTry,nconds);
for n = 1:nNBinsToTry
    nbins = nbinsToTry(n);
%     [~,binEdges] = labelByPercentile(input(OSInds),nbins);
	[~,binEdges] = labelByPercentile(balancedInput,nbins);
    [~,~,inputLabels] = histcounts(input,binEdges);
    inputLabel_bySucc_byCond = groupDataByLabel(inputLabels,[successLabels condLabels]);
    inputLabelPMFs_bySucc_byCond = cellfun(@(x) ...
        histcounts(x,0.5:nbins+0.5,'normalization','probability'),...
        inputLabel_bySucc_byCond,'uniformoutput',false);
    sKLD_byNBins_byCond(n,:) = cellfun(@(x,y) sKLD(x,y),...
        inputLabelPMFs_bySucc_byCond(1,:),inputLabelPMFs_bySucc_byCond(2,:));
end; clear n
sKLD_byCond = median(sKLD_byNBins_byCond);
end

