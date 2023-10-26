function [output] = removeMeanFRByDay(input,lastBins)
%%% This function takes in the input matrix of firing rates (assumed to
%%% comprise of all days with NaNs present when units are not in specific
%%% days) and the given last trial for each day and returns back a firing
%%% rate matrix with all of the unit FR values after subtracting their
%%% mean.
%%%
%%% inputs:
%%% input - a matrix of FRs (with NaNs at unit locations with no FR for 
%%% that day) that is size trials x units
%%% lastTrials - vector with indices of the last trial for each day
%%%
%%% outputs:
%%% output - matrix of FRs, same size as input, but now w/ mean over the
%%% day removed from each unit (like z-scoring but w/o division by std)
%%% Adam Smoulder, 10/10/18 (edited 12/14/18)
%%%

output = nan(size(input));
lastBins = [0 lastBins]; % just makes the loop easier to write
for i = 1:(length(lastBins)-1)
    indicesForDay = (lastBins(i)+1):lastBins(i+1);
    output(indicesForDay,:) = ...
            input(indicesForDay,:) - mean(input(indicesForDay,:));
end; clear i;

end

