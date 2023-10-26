function [output] = zscoreFRByDay(input,lastTrials)
%%% This function takes in the input matrix of firing rates (assumed to
%%% comprise of all days with NaNs present when units are not in specific
%%% days) and the given last trial for each day and returns back a firing
%%% rate matrix with all of the unit FR values z-scored within their given
%%% day.
%%%
%%% inputs:
%%% input - a matrix of FRs (with NaNs at unit locations with no FR for 
%%% that day) that is size trials x units
%%% lastTrials - vector with indices of the last trial for each day
%%%
%%% outputs:
%%% output - matrix of FRs, same size as input, but now z-scored
%%% Adam Smoulder, 10/10/18
%%%

output = nan(size(input));
lastTrials = [0 lastTrials]; % just makes the loop easier to write
for i = 1:(length(lastTrials)-1)
    indicesForDay = (lastTrials(i)+1):lastTrials(i+1);
    output(indicesForDay,:) = zscore(input(indicesForDay,:));
end; clear i;

end

