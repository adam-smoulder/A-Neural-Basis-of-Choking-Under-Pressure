function [y,w] = twoPointInterpolation(x,x1,x2,y1,y2)
% This function uses a weighted average to interpolate the value "y"
% associated with input "x" and its two nearest points, (x1,y1) and
% (x2,y2). This average is weighted inversely proportional to the distance
% of x from x1 and x2. This is a 2-point linear interpolation. Sure you
% could do splines or something, but often this suffices. Assumes x-values
% are all scalar (y's can be whatever so long as they're the same size).
%
% Inputs:
% - x: scalar of the input for which we want an interpolated output value
% - x1: scalar of input for nearest data point below x
% - x2: scalar of input for nearest data point above x
% - y1: output corresponding to input x1
% - y2: output corresponding to input x2
%
% Outputs:
% - y: interpolated output corresponding to x
% - w: [2 x 1] vector of weights used on y1 and y2
%
% This was made mostly for the fact that we have position (y values) 
% sampled only at certain timepoints (x values), and we want to estimate a
% position for a time between our samples.
%
% Adam Smoulder, 3/25/20

dists = abs([x1 ; x2] - x); % distance from x
w = flip(dists)/sum(dists); % weights are inversely proportional to distance
y = w(1)*y1 + w(2)*y2;      % I write it this way incase y has weird sizes

end

