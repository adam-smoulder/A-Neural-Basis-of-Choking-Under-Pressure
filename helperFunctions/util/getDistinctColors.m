function [colors] = getDistinctColors(varargin)
%%% getDistinctColors simply returns back a 1x15 cell vector of 1x3 double
%%% arrays, containing the RGB triplets of a bunch of colors that are pretty
%%% distinct from each other
%%%
%%% inputs:  none needed, but use name value pairs for the following:
%%% - 'SELECT_ORDER' - integer for a 'preset' order; defined below
%%% - 'RANDOM_ORDER' - default false; if true, randomizes color order
%%% - 'NUM_COLORS' - integer to only return back however many colors
%%%                  instead of all
%%%
%%% outputs:
%%% - colors - 1x14 cell vector of 1x3 arrays, containing the RGB triplets of a bunch of colors
%%% Adam Smoulder, 9/19/18

% take input arguments
SELECT_ORDER = 0;   % 0 = default preset
RANDOM_ORDER = 0;   % 0 = not random
NUM_COLORS = 0;     % 0 reads as all colors, just incase we add more!
warnOpts(assignOpts(varargin));


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
if RANDOM_ORDER
    colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
    colors = colors(randperm(length(colors)));
elseif SELECT_ORDER
    switch SELECT_ORDER
        case 0 % my default - R G B M O I C Y T B... also what I use for reward sizes
            colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
        case 1 % rough rainbow
            colors = [red orange yellow green cyan teal blue ...
                indigo navy plum maroon brown gray black];
        case 2 % reverse rainbow
            colors = [ red orange yellow green cyan blue teal...
                indigo navy plum maroon brown gray black];
            colors = flip(colors);
        case 3 % same as default but starts with black
            colors = [black red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray];
        case 4 % for tubes; two similar colors beside each other
            colors = [red orange green lightgreen blue cyan plum magenta ...
              maroon brown  gray black];
        case 5 % easy rainbow that starts with black with teal b/w green and blue (9 colors total)
            colors = [black red orange yellow green teal blue indigo plum];
        case 6 % used for punishments
            colors = [yellow orange cyan indigo];
        case 7 % for the 9 failure types: failStatusNames = {'MS','FS','Cheat','NA','OS','US','Slow','Drift','Jitter'};
            colors = {[0.65 0.65 0.65],[0 0.7 0.7],[0 1 1],0.9*[0.1 0.4 0.1],[0.1 0.55 0.1],[0.1 0.77 0.35],[0 1 0],[0.6 0 0.6],[1 0 1]};
        case 8 % NICK'S REWARD COLORS from the behavior paper; also adds 4 more after
            colors = [{[1,0,0],[1,0.647,0],[0,0,1],[0,0,0]},yellow,green,cyan,magenta];
        case 9 % 20 distinct colors, taken from some random website
            colors = [{[0, 0, 0]/255},{[230, 25, 75]/255}, {[60, 180, 75]/255}, {[255, 225, 25]/255}, {[0, 130, 200]/255}, {[245, 130, 48]/255}, {[145, 30, 180]/255}, {[70, 240, 240]/255}, {[240, 50, 230]/255}, {[210, 245, 60]/255}, {[250, 190, 190]/255}, {[0, 128, 128]/255}, {[230, 190, 255]/255}, {[170, 110, 40]/255}, {[255, 250, 200]/255}, {[128, 0, 0]/255}, {[170, 255, 195]/255}, {[128, 128, 0]/255}, {[255, 215, 180]/255}, {[0, 0, 128]/255}, {[128, 128, 128]/255}, {[50, 130, 240]/255}];
        case 10 % a nice 8 color wheel with constant blue
            disp('ADAM YOU SHOULDNT BE USING THIS COLOR WHEEL FOR DIRECTIONS ANYMORE 2/10/23 use case 16')
            b = 0.5;
            colors = {[1 .5 b],[.75 .75 b],[.5 1 b],[.25 .75 b],[0 .5 b],[0.25 0.25 b],[0.5 0 b],[0.75 0.25 b]};
        case 11 % another nice 8 color wheel with constant green
            g = 0.38;
            colors = {[1 g .5],[.75 g .75 ],[.5 g 1 ],[.25 g .75],[0 g .5 ],[0.25 g 0.25],[0.5 g 0 ],[0.75 g 0.25]};
        case 12 % for the rare-large 5 reward control
            colors = [{[1,0,0],[1,0.647,0],[0,0,1],0.75*cyan{1},[0,0,0]},yellow,green,magenta];
        case 13 % for the common-jackpot 5 reward control
            colors = [{[1,0,0],[1,0.647,0],[0,0,1],0.75*magenta{1},[0,0,0]},yellow,green,cyan];
        case 14 % 4 rewards but gentler
            colors = {[1,0.3,0.3],[1,0.747,0.3],[0.3,0.3,1],[0.3,0.3,0.3]};
        case 15 % green shading gradient, 8 colors
            colors = {[183 255 191]/255,[149 249 133]/255,[77 237 48]/255,[38 215 1]/255,[0 190 1]/255,[0 161 8]/255,[0 125 10]/255,[0 80 20]/255};
        case 16 % same as case 10 but reordered
            b = 0.5;
            colors = {[1 .5 b],[0.75 0.25 b],[0.5 0 b],[0.25 0.25 b],[0 .5 b],[.25 .75 b],[.5 1 b],[.75 .75 b]};
        otherwise
            colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
            colors = colors(randperm(length(colors)));
    end
else
    colors = [red green blue magenta orange indigo cyan yellow teal brown...
              maroon navy gray black];
end

% assign proper number of colors, then output!
numcolors = ternop(NUM_COLORS==0, length(colors), NUM_COLORS); 
colors = colors(1:numcolors);
end

