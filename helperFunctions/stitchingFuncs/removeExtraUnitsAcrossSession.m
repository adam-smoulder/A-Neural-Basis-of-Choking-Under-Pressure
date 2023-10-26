function [data] = removeExtraUnitsAcrossSession(data)
%%% Across some sessions, a few units drop and come back; we want to
%%% conserve only units that are present across the whole session. This
%%% prunes trails for a given dataset such that this is the case!
%%%
%%% 'data' is assumed to be a N-length structure (for N trials) with a field
%%% called 'TrialData' - 'TrialData', in turn, should have a field
%%% 'spikes', which is length M (M number of units = channels & sorts), and
%%% has a field called 'sort'
%%%
%%% Adam Smoulder, 10/3/18



% get the conserved units (I call it "minimum set")
ucIds = cell(1,length(data));  % unit-channel ids
for i = 1:length(data)
    % first, get the "ids" for the trial, meaning their channels and sorts
    currentIds = 10*[data(i).TrialData.spikes.channel]...
        +[data(i).TrialData.spikes.sort]; % unique bc only 96 chans
    ucIds{i} = currentIds;
    
    % next, intersect this with previous trials (to get minimum set)
    if i==1
        minimumSet = currentIds;
    else
        minimumSet = intersect(minimumSet, currentIds);
    end
end; clear i;


% now apply the minimum set to each
for i = 1:length(data)
    currentIds = ucIds{i};
    
    % find what differences are present between minimum set and current
    superfluousIds = setdiff(currentIds, minimumSet);
    
    % find and remove these ids in the data
    if ~isempty(superfluousIds)
        for j = length(superfluousIds):-1:1 % backwards to preserve trial order
            indexToRemove = find(currentIds == superfluousIds(j));
            data(i).TrialData.spikes(indexToRemove) = [];
            % note, this does NOT update waveforms... 
            % we may have to do that in the future
        end; clear j;
    end
end; clear i;


end

