function [dist] = eucDist(varargin)
% Finds the Euclidean distance between two points/vectors, X1 and X2.
% Assumes they're of the same dimensionality. If the inputs are matrices of
% the same dimensionality (N x M), the output is a row vector that finds the
% distances between the points along the first dimension (1 x M output)
%
% If you pass in only 1 matrix, the distance will be sequentially found
% along each dimension. THIS ONLY CURRENTLY WORKS FOR 2D MATRICES!
%
% Really not sure why MATLAB doesn't have this built in somewhere...
%
% Adam Smoulder, 6/3/19
switch nargin
    case 1
        X = varargin{1};
        xSize = size(X);
        dist = nan(xSize(1)-1,1);
        for i = 1:(xSize(1)-1) % we're assuming X is <= 7 dimensional
            dist(i) = sqrt(sum((X(i,:,:,:,:,:,:)-X(i+1,:,:,:,:,:,:)).^2));
        end; clear i
    otherwise
        X1 = varargin{1};
        X2 = varargin{2};
        dist = sqrt(sum((X1-X2).^2));
end

end

