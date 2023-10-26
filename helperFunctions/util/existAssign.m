function [structureString] = existAssign(structureString,variableString)
% if a variable by the name of variableString exists in the workspace, it
% will be assigned to structure of name structureString. Otherwise, nothing 
% will. The structure is then returned. 
%
% Adam Smoulder, 2/18/19

if ismember(variableString,evalin('base','who'))
    evalString = [structureString '.' variableString ' = ' variableString ';'];
    evalin('base',evalString)
end

end

