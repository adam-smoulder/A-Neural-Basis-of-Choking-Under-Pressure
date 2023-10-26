function [output] = shannonEntropy(data)
% Calculates Shannon entropy of the given (joint) probability matrix
%
% Inputs:
% - data: Either the (joint) probability distribution of the data or some
% potential distribution or counts or anything really - for the purposes of
% entropy, we treat all data as equal and vectorize + normalize it such
% that it becomes a big (joint) probability distribution. For your
% organization, I'd recommend you still make it a probability distribution
% before input.
%
% Outputs:
% - output: The Shannon entropy of the data's (joint) probability 
% distribution
%
% Adam Smoulder, 4/9/19

% normalize if not already done
data = data./sum(data(:)); 

% vectorize it
data = data(:);

% find entropy
output = -sum(data.*log2(data+realmin)); % + realmin avoids log of 0

end

