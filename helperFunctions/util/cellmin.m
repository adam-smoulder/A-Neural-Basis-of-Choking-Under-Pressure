function [output] = cellmin(input)
% Minimum value for each cell, using cellfun. Just a bit easier to write.
%
% Adam Smoulder, 5/21/20
output = cellfun(@(x) min(x), input,'uniformoutput',false);
end

