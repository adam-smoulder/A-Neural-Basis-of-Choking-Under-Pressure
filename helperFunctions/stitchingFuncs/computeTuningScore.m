function [tuningScore] = computeTuningScore(tuning)
% An extension to the Fraser/Schwartz code pack that takes into account
% tuning via a "tuning score". The score formula to compare the similarity
% in two tuning curves is:
%    sum((abs(mu1 - mu2)/sqrt(var1+var2))
% where the sum is over conditions (e.g. tube*epoch)
%
% Inputs:
% - tuning: ndays x 1 cell structure, where each cell contains a cell array
% that's nunits x 1 for the given day, where each cell is nconds x 2. The
% first value of the 2 is the mean for the condition, the second is the
% variance.
%
% Outputs:
% - tuningScore: cell with a nunits1 x nunits2 double array containing the
% tuning similarity for each pair of units
%
% ONLY CURRENTLY DESIGNED FOR TWO DAYS AT A TIME
% Adam Smoulder, 6/3/19

% get stuff we'll need
tuning1 = tuning{1};
tuning2 = tuning{2};
tuningScore = zeros(length(tuning1),length(tuning2));

% iterate over pairwise comparisons and get that unit score
for i = 1:length(tuning1)
    mean1 = tuning1{i}(:,1);
    var1 = tuning1{i}(:,2);
    for j = 1:length(tuning2)
        mean2 = tuning2{j}(:,1);
        var2 = tuning2{j}(:,2);
        tuningScore(i,j) = sum(abs(mean1-mean2)./sqrt(var1+var2));
    end; clear j
end; clear i

tuningScore = {tuningScore};

end

