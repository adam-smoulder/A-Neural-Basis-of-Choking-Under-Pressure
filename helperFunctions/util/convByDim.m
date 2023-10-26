function [convMatrix] = convByDim(matrix,kernel,dim,varargin)
% For some great reason, you can't just take a 1D kernel and convolve it
% along each row/column of matrix in MATLAB's conv function. Here, we
% address that.
%
% Inputs:
% - matrix: [? ... x n1 x ... ?] matrix that we'll be convolving along. The
% dimension specified by "dim" is the one we'll convolve along (n1 in
% length)
% - kernel: [n2 x 1] vector with the convolution kernel
% - dim: integer, the dimension that we're going to convolve along for
% matrix
%
% Optional Inputs:
% - convArg: (default 'full') additional SHAPE argument to use with conv;
% other options are 'same' and 'valid' (use "help conv" for details)
%
% Outputs:
% - convMatrix: matrix, but now with one dimension convolved with kernel.
% Size may be larger than input based on convArg
%
% Adam Smoulder, 4/24/20

convArg = 'full';
if length(varargin)==1
    convArg = varargin{1};
end

% First, rearrange matrix such that dim is at the front, then flatten it
origSize = size(matrix);
ndims = length(origSize);
permOrder = [dim setdiff(1:ndims,dim)];
permMatrix = permute(matrix,permOrder);
permSize = size(permMatrix);
flatMatrix = permMatrix(:,:);

% Now do the convolutions! Different convArgs will yield different sized
% outputs, so we'll adapt based on the first one
for i = 1:size(flatMatrix,2)
    curVector = flatMatrix(:,i);
    if i == 1 % we need to setup the full matrix
        exRow = conv(curVector,kernel,convArg);
        convFlatMatrix = nan([length(exRow),size(flatMatrix,2)]);
    end
    
    % Now we actually do the stuff!
    curVector(isnan(curVector)) = [];                          % skip nans
    curConvVector = conv(curVector,kernel,convArg);            % do conv
    convFlatMatrix(1:length(curConvVector),i) = curConvVector; % assign
end; clear i


% Reshape and repermute the convFlatMatrix to get it back to normal shape
% First, set up sizes of matrices 
convPermSize = [size(convFlatMatrix,1) permSize(2:end)];
convSize = origSize;
convSize(dim) = size(convFlatMatrix,1);

% Next, reshape to get it back to the permuted shape
convPermMatrix = reshape(convFlatMatrix,convPermSize);

% and finally, permute it back to being the original shape
depermOrder = [2:dim 1 (dim+1):ndims]; % weird looking, but correct
convMatrix = permute(convPermMatrix,depermOrder);

end

