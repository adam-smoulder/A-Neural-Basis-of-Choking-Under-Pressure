function plotCircle_adam(varargin)
% Plots a circle. Yeehaw. Default is radius = 1, mean = [0 0], color is
% black, linewidth is 1. Only use for 2D plots
%
% Inputs:
% - radius: radius of the circle (scalar)
% - center:  where the center of the circle is, [2x1] or [1x2]
% - color:  the color of the circle ([1x3] or [3x1] for RGB, on (0,1))
% - linewidth:  thickness of the circle
%
% Output: None.
%
% Assumes that a figure is already made and already has "hold on"
%
% Adam Smoulder, 3/23/2020

% defaults
radius = 1;
center = [0 0];
color = [0 0 0];
linewidth = 1;

% inputs
if length(varargin) >= 1
    radius = varargin{1};
end
if length(varargin) >= 2
    center = varargin{2};
end
if length(varargin) >= 3
    color = varargin{3};
end
if length(varargin) >= 4
    linewidth = varargin{4};
end

% plot the damn circle
npts = 100;
theta = linspace(0,2*pi,npts);
xvals = radius*cos(theta)+center(1);
yvals = radius*sin(theta)+center(2);
plot(xvals,yvals,'-','color',color,'linewidth',linewidth)

end

