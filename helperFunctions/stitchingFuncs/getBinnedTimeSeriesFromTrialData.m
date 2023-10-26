function [binnedData, binTime] = getBinnedTimeSeriesFromTrialData(trialData,binSize,startEvent,startDelay,endEvent,endDelay)
% This function takes in trialData (a sructure output after stitching) and
% gets a time-series of non-overlapping binned factor scores for each trial
% based on the input specifications of bin size/timing.
%
% Inputs:
% - trialData: [ntrials x 1] struct with all data and timing labels
% - binSize: scalar (ms) of what size the bin should be
% - startEvent: string stating which event to align the start of bin 1 to.
% See getSingleDataBinFromTrialData for options of events that you can
% align the time series to. If the event time is "nan" for a
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
%
% Outputs:
% - binnedData: [ntrials x 1] cell array with the binned
% counts/trajectories/whatever for each trial.
% - binTime: [ntime x 1] vector of the time vector for the last trial that was
% binned. If you're binning the same interval every time, this provides an
% easy output to take advantage of.
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
% Another note: I've only ever really tried this with 1ms and then did a
% smoothing kernel afterwards. I think I tested it when I wrote it to check
% if other bin sizes work? 
%
% Adam Smoulder, 1/14/21 (edit 9/30/21)

ntrials = length(trialData);
binnedData = cell(ntrials,1);
binTime = [];
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
    
    % Bin the data and store it
    time = trialData(i).time;
    trialMatrix = [];
    for j = 1:nbins
        startInd = find(time >= binEdges(j),1);
        endInd = find(time >= binEdges(j+1),1)-1;
        trialMatrix(end+1,:) = mean(trialData(i).neuralData.factorScores(:,startInd:endInd),2)';
    end; clear j
    
    % If it's the first trial, do the binned time as well
    if isempty(binTime)
        for j = 1:nbins
            startInd = find(time >= binEdges(j),1);
            endInd = find(time >= binEdges(j+1),1)-1;
            binTime(end+1) = mean([time(startInd) time(endInd)]);
        end; clear j
        binTime = binTime-startEventTime;
    end
    
    binnedData{i} = trialMatrix;
end; clear i

end

