function [unitCid] = cidFromId(unitId)
%  Calculates unit conmbined ID # (cid) from the unit id. Shouldn't need
%  changed, unless it's for efficiency or perhaps larger arrays?
%
%  inputs:
%  - unitId: n x 3 matrix, where for n units, the first column is the day,
%  second column is the channel, and third column is the sort
%
%  outputs:
%  - unitCid: n x 1 vector of integers, calculated with some simple cipher
%  from the unitId
%
% Adam Smoulder, 9/18/18 (edited: 2/14/19)

if isempty(unitId)
    unitCid = [];
else
    % ID from CID is (just for reference):
    %unitId = [floor(unitCid/1000) mod(floor(unitCid/10),floor(unitCid/1000)*100) mod(unitCid,floor(unitCid/10)*10)];
    
    %unitCid = unitId*[1000 10 1]';
    unitCid = sum(bsxfun(@times, unitId, [1000 10 1]),2); % using this for matlab 2015 compatability
end

end

