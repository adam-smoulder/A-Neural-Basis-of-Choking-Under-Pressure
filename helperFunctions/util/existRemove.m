% This function takes in a structure (myStruct) and a cell array of
% putative fields in the structure that you want to remove (myFields). It
% then runs through each of the fields and removes them if they exist in
% the structure.
%
% Inputs:
% - myStruct: A structure that you want to remove fields from
% - myFields: 1D array/vector of cells, where each cell contains a string
% with the name of a field you want to remove
%
% Optional Input:
% - displayText: Display if the field was not present on the structure
% inputted. (Default: true)
%
% Outputs:
% - myStruct: Your structure, but now with myFields removed
%
% I made this because if you just use rmfield, it throws an error if the
% field doesn't exist (why it does this instead of a warning is beyond me).
%
% Adam Smoulder, 1/22/20 (edit 9/7/20)

function [myStruct] = existRemove(myStruct,myFields,varargin)

if length(varargin) ~=0
    displayText = varargin{1};
else
    displayText = true;
end

for i = 1:length(myFields) % for each field
    if isfield(myStruct,myFields{i}) % if this field is in the structure
        myStruct = rmfield(myStruct,myFields{i}); % remove it
    else % otherwise, say so
        if displayText
            disp(['No field in structure called:  ' myFields{i}])
        end
    end
end; clear i
end

