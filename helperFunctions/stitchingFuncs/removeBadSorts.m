function [data] = removeBadSorts(data)
%%% Some data has a bad sort in it for some reason... this function removes
%%% them. They should all be sort id 31, and so we just kill em. NOTE: this
%%% does not resave data; this is just meant to be run to sift them out
%%% before resaving in another context.
%%%
%%% inputs (and outputs in this case):
%%% 'data' is assumed to be a N-length structure (for N trials) with a field
%%% called 'TrialData' - 'TrialData', in turn, should have a field
%%% 'spikes', which is length M (M number of units = channels & sorts), and
%%% has a field called 'sort'
%%%
%%% Adam Smoulder, 9/1/18 (edit 10/3/18)

badTrials = [];
for trialCount = 1:length(data)
    if ~isempty(data(trialCount).TrialData.spikes)
        keepers = find([data(trialCount).TrialData.spikes.sort]~=31);
        data(trialCount).TrialData.spikes = data(trialCount).TrialData.spikes(keepers);
    else
        badTrials = [badTrials trialCount];
    end
end

data(badTrials) = [];
end

