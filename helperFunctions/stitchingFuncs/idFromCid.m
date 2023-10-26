function [unitId] = idFromCid(unitCid)
%  Calculates unit conmbined ID # (cid) from the unit id. Shouldn't need
%  changed, unless it's for efficiency or perhaps larger arrays?
%
%  inputs:
%  - unitCid: n x 1 vector of integers, calculated with some simple cipher
%  from the unitId
%
%  outputs:
%  - unitId: n x 3 matrix, where for n units, the first column is the day,
%  second column is the channel, and third column is the sort
%
%  Adam Smoulder, 2/14/19


if isempty(unitCid)
    unitId = [];
else
    %unitCid = unitId*[1000 10 1]'; cID to ID
    unitId = [floor(unitCid/1000) mod(floor(unitCid/10),floor(unitCid/1000)*100) mod(unitCid,floor(unitCid/10)*10)];
end

end

