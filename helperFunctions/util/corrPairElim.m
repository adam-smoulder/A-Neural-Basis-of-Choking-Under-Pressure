function [elimUnits] = corrPairElim(spikeCorr, corrThresh)
% corrPairReduce is made to take the square correlation matrix of channels'
% spike data and output which units to remove from further analysis in
% order to not double-count. Having near identical unit firing may happen
% due to electrical shorts or things of that nature.
%
% inputs:
% - spikeCorr:  square correlation matrix, where max value is 1.0 (nunits x
% nunits)
% - corrThresh:  threshold corr value to call two units "same" or not
% (recommended 0.8ish, though lower is more conservative). 
% 
% outputs:
% - elimUnits:  row vector of indices indicating which units to elimate
%
% Adam Smoulder, 1/7/19, edit 11/11/19
%

nunits = size(spikeCorr,1); % number of units

% find where correlation is above threshold (aside from self-correlation)
nonSelfCorr = spikeCorr-diag(diag(spikeCorr));
currBadInds = find(nonSelfCorr > corrThresh);

% arrange it in a "pairwise" format
badPairs = [ceil((currBadInds-1)/nunits) mod((currBadInds-1),nunits)+1];
elimUnits = [];
while ~isempty(badPairs)
    unitToElim = mode(badPairs(:));  % find unit w/ most high correlations
    elimUnits = [elimUnits unitToElim]; % add it to the list to kill later
    elimRows = mod(find(badPairs == unitToElim)-1,size(badPairs,1))+1;
    badPairs(elimRows,:) = []; % remove it from the current list
end

% additionally, if self-corr is nan, it's gonna be a bad unit (not firing
% enough consistently) % actually don't default use this, turns out it's
% inconsistent and will change between runs...let other methods remove
% these bad units
% nanBads = find(isnan(diag(spikeCorr)));
% elimUnits = [elimUnits nanBads'];

elimUnits = sort(elimUnits); % sort for good measure
end

