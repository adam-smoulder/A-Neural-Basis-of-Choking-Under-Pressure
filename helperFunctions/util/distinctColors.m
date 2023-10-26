function [colors] = distinctColors(varargin)
% getDistinctColors simply returns back a 1x15 cell vector of 1x3 double
% arrays, containing the RGB triplets of a bunch of colors that are pretty
% distinct from each other
%
% inputs:  none needed, but you can input the following in order:
% - option:  scalar from 0-7. One of the preset color options. Default=0
% - ncolorout:  scalar of how many colors are desired out (max of 13 for 
% most orders). Default=all
% - randomize;  boolean of if the color order should be randomized.
% Supercedes the "option" parameter.
%
% outputs:
% - colors - 1x13 cell vector of 1x3 arrays, containing the RGB triplets of
% a bunch of colors
%
% Adam Smoulder, 9/19/18

% take input arguments
if nargin == 3
    randomize = varargin{3};
    ncolorout = varargin{2};
    option = varargin{1};
elseif nargin == 2
    randomize = 0;
    ncolorout = varargin{2};
    option = varargin{1};
elseif nargin == 1
    randomize = 0;
    ncolorout = 0; 
    option = varargin{1};
else % nargin == 0 or some odd input with more
    randomize = 0; % not a random order!
    ncolorout = 0; % just does all colors
    option = 0; % default preset
end


% define color rgb values:
red =       {[230, 25, 75]/255};
green =     {[60, 180, 75]/255};
lightgreen= {[60, 180, 75]/190};
blue =      {[0, 130, 200]/255};
plum =      {[240, 50, 230]/(255*1.5)};
magenta =   {[240, 50, 230]/255};
yellow =    {[255, 225, 25]/255};
orange =    {[245, 130, 48]/255};
indigo =    {[75,0,130]/255};
cyan =      {[70, 240, 240]/255};
teal =      {[0, 128, 128]/255};
brown =     {[170, 110, 40]/255};
maroon =    {[128, 0, 0]/255};
navy =      {[0, 0, 128]/255};
gray =      {[47, 79, 79]/255};
black =     {[0,0,0]/255};


% change the order if desired
if randomize
    colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
    colors = colors(randperm(length(colors)));
elseif option
    switch option
        case 0 % my default - R G B M O I C Y T B...
            colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
        case 1 % rough rainbow
            colors = [red orange yellow green teal blue ...
                indigo navy plum maroon brown gray black];
        case 2 % reverse rainbow
            colors = [ red orange yellow green cyan blue teal...
                indigo navy plum maroon brown gray black];
            colors = flip(colors);
        case 3 % same as default but starts with black
            colors = [black red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray];
        case 4 % for pairs; two similar colors beside each other
            colors = [red orange green lightgreen blue cyan plum magenta ...
              maroon brown  gray black];
        case 5 % easy rainbow that starts with black with teal b/w green and blue (9 colors total)
            colors = [black red orange yellow green teal blue indigo plum];
        case 6 % some other quick distinct 4
            colors = [yellow orange cyan indigo];
        case 7 % 20 distinct colors, taken from some random website
            colors = [{[0, 0, 0]/255},{[230, 25, 75]/255}, {[60, 180, 75]/255}, {[255, 225, 25]/255}, {[0, 130, 200]/255}, {[245, 130, 48]/255}, {[145, 30, 180]/255}, {[70, 240, 240]/255}, {[240, 50, 230]/255}, {[210, 245, 60]/255}, {[250, 190, 190]/255}, {[0, 128, 128]/255}, {[230, 190, 255]/255}, {[170, 110, 40]/255}, {[255, 250, 200]/255}, {[128, 0, 0]/255}, {[170, 255, 195]/255}, {[128, 128, 0]/255}, {[255, 215, 180]/255}, {[0, 0, 128]/255}, {[128, 128, 128]/255}, {[50, 130, 240]/255}];
        otherwise
            colors = colors(randperm(length(colors)));
    end
else
    colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
end

% assign proper number of colors, then output!
if ncolorout==0, ncolorout = length(colors); end
colors = colors(1:ncolorout);
end

