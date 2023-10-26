function [infoVec] = optSubplotParams(n,flip)
% optSubplot returns the best values for splitting up your subplots and
% stuff based on how many panels there should be total (n). The split is
% "optimized" first based on not leaving blanks, then on squareness, and
% and also on accounting for just raw # of panels. Works on the assumption
% that you have more room across than vertical, but if you don't wanna do
% that, use flip.
% 
% inputs:  n - number of panels to show, flip - 0 if normal, 1 if switch
% the numbers of the output rows and cols
%
% outputs: infoVec - contains different info at each index:
%   1 - number of rows
%   2 - number of columns
%   3 - where the x-axis label should be
%   4 - where the y-axis label should be
%   5 - where the title should be
%   6 - where a legend should be (upper right panel)
%   
%   Adam Smoulder, 12/17/18

% I originally tried figuring out the math to make this automated, but it
% honestly was such a pain that I just pretty much made a dictionary here,
% and if it's above 40 then it picks bs-y
infoVec = zeros(1,5);

switch n
    case 1
        infoVec(1:2) = [1 1];
    case 2
        infoVec(1:2) = [1 2];
    case 3
        infoVec(1:2) = [1 3];
    case 4
        infoVec(1:2) = [2 2];
    case 5
        infoVec(1:2) = [1 5];
    case 6
        infoVec(1:2) = [2 3];
    case 7
        infoVec(1:2) = [2 4];
    case 8
        infoVec(1:2) = [2 4];
    case 9
        infoVec(1:2) = [3 3];
    case 10
        infoVec(1:2) = [2 5];
    case 11
        infoVec(1:2) = [3 4];
    case 12
        infoVec(1:2) = [3 4];
    case 13
        infoVec(1:2) = [3 5];
    case 14
        infoVec(1:2) = [3 5];
    case 15
        infoVec(1:2) = [3 5];
    case 16
        infoVec(1:2) = [4 4];
    case 17
        infoVec(1:2) = [3 6];
    case 18
        infoVec(1:2) = [3 6];
    case 19   
        infoVec(1:2) = [4 5];
    case 20
        infoVec(1:2) = [4 5];
    case 21
        infoVec(1:2) = [4 6];
    case 22
        infoVec(1:2) = [4 6];
    case 23
        infoVec(1:2) = [4 6];
    case 24
        infoVec(1:2) = [4 6];
    case 25
        infoVec(1:2) = [5 5];
    case 26
        infoVec(1:2) = [4 7];
    case 27
        infoVec(1:2) = [4 7];
    case 28
        infoVec(1:2) = [4 7];
    case 29
        infoVec(1:2) = [5 6]; 
    case 30
        infoVec(1:2) = [5 6];
    case 31
        infoVec(1:2) = [5 7];
    case 32
        infoVec(1:2) = [5 7];
    case 33
        infoVec(1:2) = [5 7];
    case 34
        infoVec(1:2) = [5 7];
    case 35
        infoVec(1:2) = [5 7];
    case 36
        infoVec(1:2) = [6 6];
    case 37
        infoVec(1:2) = [5 8];
    case 38
        infoVec(1:2) = [5 8];
    case 39
        infoVec(1:2) = [5 8];
    case 40
        infoVec(1:2) = [5 8];
    otherwise
        infoVec(1:2) = [ceil(n/8) 8];
end

% switch em if desired
if flip
    infoVec(1:2) = infoVec([2,1]);
end

% ideal axis label locations
infoVec(3) = floor(infoVec(2)/2)+infoVec(2)*(infoVec(1)-1)+1; % x-axis
infoVec(4) = infoVec(2)*floor((infoVec(1)-1)/2)+1; % y-axis
infoVec(5) = floor(infoVec(2)/2)+1; % title
infoVec(6) = infoVec(2);