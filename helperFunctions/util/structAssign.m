function [structure] = structAssign(structure,fields,varargin)
% This function takes a "structure" with N entries and assigns each entry
% with the values in "fields" (be them new fields or existing ones). Each
% field in fields is assumed to also have N entries.
%
% Inputs:
% - structure: the structure you want to assign the field values to, where
% length(structure) = N
% - fields: a cell array that is nfields*2 in length. The first cell
% should contain the string with the name of the field you wish to assign.
% The second cell should have N rows (as its first dimension length) and
% contain the values you want to be assigned to the field for each entry in
% structure. Then the third cell is your second field you want to assign,
% the fourth cell contains its values, and so on.
%
% Optional Input:
% - unboxCells: boolean (default true); if a field's values are cells (e.g.
% a N x 1 cell array), should we open up the cell for each index and assign
% its contents to the structure? If you want to keep the fields in the
% structure as cells, say "false"
%
% Output:
% - structure: same as input structure, but with updated fields
%
%
% Ex: If I have a 300 entry structure (myStructure) and I wanted to add 
% some weird fields to it:
%
% myFields = {'randomIntField',randi(10,[300 1]),...
%             'randomNumMatrix',randn([300 10]),...
%             'randomCells',cell(100,1)};
% updatedStructure = structAssign(myStructure, myFields, true);
%
% Adam Smoulder, 3/24/2020

% parse optional arguments
unboxCells = true;
if length(varargin) > 1
    unboxCells = varargin{1};
end

% make sure the structure is either the same length as the first field or
% it is empty
assert(isempty(structure) || length(structure)==size(fields{2},1))

% assign the fields
for i = 1:length(fields)/2
    curFieldName = fields{2*i-1};
    curFieldVals = fields{2*i};
    for j = 1:length(curFieldVals)
        structure(j).(curFieldName) = curFieldVals(j,:,:,:,:,:); % assumes your data is, at max, a 6-dim tensor...
    end; clear j
end; clear i

end

