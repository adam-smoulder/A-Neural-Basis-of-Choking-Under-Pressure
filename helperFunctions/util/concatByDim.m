function [output] = concatByDim(matrix1, matrix2, dimToConcat)
% So you have a >2D matrix and you want to concatenate it along a given
% dimension that's not the first (you'd just [matrix1 ; matrix2]) or the
% second (use [matrix1 matrix2]). What do you do? You use concatByDim!
% Simply enter your two matrices along with the dimension you want to
% concatenate them along, and voila, they shall be concatenated.
%
% Only works if all other dimensions (aside from dimToConcat) are equal in
% size!
%
% Inputs: 
% - matrix1: your first matrix that you want earlier / in front / on top
% - matrix2: the second matrix, which you want later / in back / below
% - dimToConcat: whichever dimension you want to concatenate along
%
% Output: same size as the inputs except for along the dimToConcat
%
% Ex:
% mat1 = ones([5 4 7 2]);
% mat2 = zeros(5 4 10 2]);
% mat12 = concatByDim(matrix1,matrix2,3)
%
% mat12 will be a [5 x 4 x 17 x 2] matrix with 1s and 0s.
%
% Adam Smoulder, 7/17/19 (edited 2/17/20 to work with empty arrays as matrix1)

if isempty(matrix1) % concat matrix2 to nothing; so output = matrix2
    output = matrix2;
else
    % first move the desired dimension up front and concatenate
    dims = 1:length(size(matrix1));
    noCatDims = find(dims~=dimToConcat);
    permMat1 = permute(matrix1,[dimToConcat, noCatDims]);
    permMat2 = permute(matrix2,[dimToConcat, noCatDims]);
    permOutput = [permMat1 ; permMat2];
    
    % then, repermute dimensions to put it in the original order
    outDimOrder = [noCatDims(1:dimToConcat-1)+1   1   noCatDims(dimToConcat:end)];
    output = permute(permOutput,outDimOrder);
end

end

