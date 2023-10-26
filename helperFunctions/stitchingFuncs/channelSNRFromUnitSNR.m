function [channelSNR, channelLabels] = channelSNRFromUnitSNR(unitSNR,unitCids)
% Spits out SNR by channel as opposed to unit. unitSNR should be a column
% vector with SNR for each unit in the order indicated by unitCids
%
% Adam Smoulder, 6/14/19

unitIds = idFromCid(unitCids);
channelLabels = unitIds(:,2);
channels = unique(channelLabels);
channelSNR = nan(length(channels),1);

for c = 1:length(channels)
    curInds = channelLabels==channels(c);
    channelSNR(c) = mean(unitSNR(curInds));
end; clear c

end

