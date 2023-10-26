function [output] = ternop(condition,a,b)
% Works like C# or Java ternary operator; if 'condition' is true, the
% function returns input 'a', whereas if its false, input 'b' will be
% returned.
%
% ex:  
% cond = 3<5;
% ternop(cond, 'Yes 3 is less than 5', 'Whoa 3 isn't less than 5...?')
%
% Adam Smoulder, 9/11/18 (edit 4/1/19)
%

% if(condition), output = a; else, output = b; end %old version only worked for 1 condition...

% new version expands it to as many as we want! If output, a, and b
% are of the same size as the input condition, it goes element by element
if numel(a)==numel(condition) && numel(b)==numel(condition)
    output = nan(size(condition));
    for i = 1:length(condition(:))
        if(condition(i))
            output(i) = a(i);
        else
            output(i) = b(i);
        end
    end; clear i
else % old version; just evaluate the single condition
    if condition
        output = a;
    else
        output = b;
    end
end

