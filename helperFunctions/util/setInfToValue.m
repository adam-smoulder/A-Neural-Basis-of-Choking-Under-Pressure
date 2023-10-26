function [input] = setInfToValue(input,value)
% Sometimes you have an edge case where Inf shows up and it's not
% particularly meaningfully different from some other high value, but
% because it's Inf, it breaks shit. To avoid this, we use this function to
% replace every instance of Inf in the input to whatever value you want.
%
% Inputs:
% - input: your matrix of whatever size that has infs you don't like
% - value: the value you want to set all of the infs to
% 
% Outputs:
% - input: your matrix, but now with the infs replaced
%
% This is easy to do in line, but I needed to break out this function here
% so I could call it in anonymous functions on cell arrays.
%
% Adam Smoulder, 5/28/2021

input(input==inf) = value;

end

