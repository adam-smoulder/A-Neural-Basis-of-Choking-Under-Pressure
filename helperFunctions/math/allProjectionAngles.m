function [O] = allProjectionAngles(A,B)
% this takes two matrices (A and B) assumed to be of size [D x Da] and
% [D x Db] respectively and finds the angles between each of their
% constituent column vectors in a pairwise manner. The idea is mostly to
% use this to compare projection vectors across conditions (e.g. is the
% best predictive dimension for reward orthogonal or aligned to the best
% predictive dimension for reach direction?).
%
% Inputs: 
% - A: The first matrix; columns are vectors. This
% runs is most efficient technically if A is the smaller matrix of A and B
% - B: The second matrix
%
% Outputs:
% - O: [Da x Db] matrix with the angle between each column vector of
% A and each column vector of B. For instance, O(2,1) is the angle between
% the 2nd column vector of A and the 1st column vector of B
%
% Adam Smoulder, 4/5/19
%

Da = size(A,2);
Db = size(B,2);
O = nan(Da,Db);
for i = 1:Da % for each column vector in A
    curA = A(:,i);  % get current vector
    O(i,:) = acos((curA'*B)./(norm(curA)*sqrt(sum(B.^2)))); % 1 x Db vector of angles for column i of A vs. all of B's columns
end; clear i


end

