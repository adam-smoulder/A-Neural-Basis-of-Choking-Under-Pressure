function [binnedData] = getSingleDataBinFromSpikeMatrix_adapted(spikeMat,trialData,binSize,event,delay)
% This function takes in trialData (a structure output after stitching) and
% gets one datapoint for each trial based on the input specifications of
% bin size/timing.
%
% Inputs:
% - spikeMat: [ntrials x 1] cell with the spike mats you want to use (as
% opposed to the stuff in trialData)
% - trialData: [ntrials x 1] struct with all data and timing labels
% - binSize: scalar (ms) of what size the bin should be
% - event: string stating which event to align the start of the bin to.
% Options include: t1_trialStartTime, t2_targetOnsetTime, t3_goCueTime, 
% t4_reactionTime, t4b_reactionTime_20mm, t4c_reactionTime_exit, 
% t5_peakSpeedTime, t5_1_homingStartTime, t6_reachEndTime, 
% t6b_reachEndTime_20mm, t6c_reachEndTime_entry, t6_1_homingEndTime, 
% t7_postReachStateTime, t8_trialEndTime. If the event time is "nan" for a
% given trial, the output in binnedData will be "nan" for that index.
% - delay: scalar (ms) to shift the start of the bin off of the event time.
% E.g. if you wanted target onset + 50 ms, you'd use "t2_targetOnsetTime"
% and 50 as your event and delay arguments, respectively.
%
% Outputs:
% - binnedData: [ntrials x nunits] matrix containing the binned results
%
% Adam Smoulder, 1/11/21

nunits = size(spikeMat{1},1);
binnedData = nan(length(trialData),nunits);
for i = 1:length(trialData)
    eventTime = trialData(i).(event);
    if ~isnan(eventTime)
        time = trialData(i).time;
        startInd = find(time>=eventTime+delay);
        endInd = find(time>=eventTime+delay+binSize);
        binnedData(i,:) = mean(spikeMat{i}(:,startInd:endInd),2)*1000;
    end
end; clear i
end

