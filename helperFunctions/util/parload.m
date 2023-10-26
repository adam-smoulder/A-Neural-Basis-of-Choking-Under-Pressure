function [output] = parload(filename,variableName)
% load for a single variable to be used in parallel computing shenanigans
t = getCurrentTask();
output = load([filename num2str(num2str(mod(t.ID,2)+1))],variableName);
end

