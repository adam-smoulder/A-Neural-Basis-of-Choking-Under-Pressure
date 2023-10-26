function [output] = existAndTrue(theVariable)
% Exist and true checks first if the variable exists, and if so, if it's
% true.
%
% input: string of variable name ( so use existAndTrue('myVar') )
%
% output: 1 (true) if it exists and is true. 0 (false) o/w
%
% Adam Smoulder, 1/18/19

output = 0;
if evalin('base',['exist(''' theVariable ''') == 1']) % if it exists in the main workspace
    output = boolean(evalin('base',['eval(''' theVariable ''')'])); % set the output to a casted boolean
end

end

