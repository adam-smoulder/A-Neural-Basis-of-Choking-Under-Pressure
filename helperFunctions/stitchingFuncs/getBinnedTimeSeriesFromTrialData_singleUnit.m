function [binnedData] = getBinnedTimeSeriesFromTrialData_singleUnit(trialData,binSize,startEvent,startDelay,endEvent,endDelay,taskInfo)
% This function takes in trialData (a sructure output after stitching) and
% gets a time-series of non-overlapping binned spike counts for each trial
% based on the input specifications of bin size/timing.
%
% Inputs:
% - trialData: [ntrials x 1] struct with all data and timing labels
% - binSize: scalar (ms) of what size the bin should be
% - startEvent: string stating which event to align the start of bin 1 to.
% Options include: t1_trialStartTime, t2_targetOnsetTime, t3_goCueTime, 
% t4_reactionTime, t4b_reactionTime_20mm, t4c_reactionTime_exit, 
% t5_peakSpeedTime, t5_1_homingStartTime, t6_reachEndTime, 
% t6b_reachEndTime_20mm, t6c_reachEndTime_entry, t6_1_homingEndTime, 
% t7_postReachStateTime, t8_trialEndTime. If the event time is "nan" for a
% given trial, the output in binnedData will be "nan" for that index.
% - startDelay: scalar (ms) to shift the start of the first bin off of the  
% startEvent time.
% - endEvent: EITHER a string with the name of the event (similar to
% startEvent) to end the binning at, OR a scalar (ms) of how long after the
% start event to bin. This DOES TAKE INTO ACCOUNT DELAYS if you use the
% numerical option!!! Once again, if the event doesn't exist for a given
% trial, nans will be output.
% - endDelay: scalar (ms) to shift the end of the last bin off of the
% endEvent time.
% - taskInfo: structure, output from stitching stuff that we'll use to get
% the unitID mappings and number of units
%
% Outputs:
% - binnedData: [ntrials x 1] cell array with the binned
% counts/trajectories/whatever for each trial.
%
% E.g. If you wanted to get binned data for the whole delay period starting
% from -50ms before target onset until go cue with 5ms bins, you'd call: 
% getBinnedTimeSeriesFromTrialData(trialData,5,'t2_targetOnsetTime',-50,'t3_goCueTime',0);
%
% Note - if your bin sizes don't line up with the events perfectly, you'll
% just lose one bin off the end. That is, if you have 200ms bins, and your
% current bin calls for [800 1000]ms but your end event time is 950ms, the
% bin will simply not be taken.
%
%
% Adam Smoulder, 3/12/21

ntrials = length(trialData);
binnedData = cell(ntrials,1);
nunits = length(taskInfo.stitchingModel.d);
for i = 1:ntrials
    if mod(i,500)==0, disp(['Binned trial ' num2str(i)]); end
 
    % Get start time
    startEventTime = trialData(i).(startEvent);
    if isnan(startEventTime)
        binnedData{i} = nan;
        continue
    end
    startTime = startEventTime+startDelay;
    
    % Get end time
    if isnumeric(endEvent)
        endEventTime = startEventTime+endEvent;
    else
        endEventTime = trialData(i).(endEvent);
    end
    if isnan(endEventTime)
        binnedData{i} = nan;
        continue
    end
    endTime = endEventTime+endDelay;
    
    % Get bin edges
    binEdges = startTime:binSize:endTime;
    nbins = length(binEdges)-1;
    
    % Before binning, identify which units we have
    unitIDInds = taskInfo.unitIDs{trialData(i).day};
    
    % Bin the data and store it
    time = trialData(i).time;
    trialMatrix = single(nan(nbins,nunits));  % temporary; a smarter way to do this would be use int8 or something and use the max value as nan, then clear it later
    for j = 1:nbins
        startInd = find(time >= binEdges(j),1);
        endInd = find(time >= binEdges(j+1),1)-1;
        trialMatrix(j,unitIDInds) = sum(trialData(i).neuralData.spikeMatrix(:,startInd:endInd),2)'*1000/binSize; % convert to Hz
    end; clear j
    binnedData{i} = trialMatrix;
end; clear i

end

