% Assumes 2D matrix input; just finds the first 1/true for each row
% Adam Smoulder, 2/4/21

function [vals] = findFirstForEachRow(input)
    vals = nan(size(input,1),1);
    for i = 1:size(input,1)
        value = find(input(i,:),1,'first');
        if isempty(value)
            vals(i) = nan;
        else
            vals(i) = value;
        end
    end; clear i
end
