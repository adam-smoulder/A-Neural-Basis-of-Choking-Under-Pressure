function [zStats,pvals] = calculateZStatistic(correctTrials,numberTimesSeen)

% This function will calculate the z statistic from the binomial proportion 
% test for the success rates between each reward/punishment condition. It 
% uses the two-tailed binomial proportion test. The z stat needs to be 
% greater than 1.96 or less than -1.96 to have p < 0.05.
%
% Inputs:
%
% correctTrials = numRewards x nconds matrix of correct trials for each
% reward for each condition
%
% numberTimesSeen = numRewards x nconds matrix of the total number of 
% trials to each reward*cond
%
%
% Outputs:
% zStats = [nrewards x nrewards x nconds] matrix with values of [0 1 2 3]
% correpsonding to no significane, p < 0.05, p < 0.01, and p < 0.001
% 
% pvals = same size as zStats but contains the p-values for each
% comparison. This is the probability that the null hypothesis (the two
% distributions are equal) is true.
%
% Nick Pavlosky, 2018 (edited by Adam Smoulder, 4/17/19)

% critical significance value of z-stat
crit = [1.96 2.58 3.3];

% Make a matrix to store the z statistic in. 
[nrewards,nconds] = size(correctTrials);
zStats = nan(nrewards,nrewards,nconds);

% Run the binomial proportion test on the counts for each condition.
for i = 1:nrewards
    n1 = numberTimesSeen(i,:);
    p1 = correctTrials(i,:)./n1;
    
    for j = 1:nrewards
        n2 = numberTimesSeen(j,:);
        p2 = correctTrials(j,:)./n2;
        
        phat = (n1.*p1 + n2.*p2)./(n1 + n2);
        
        zStats(i,j,:) = (p1 - p2)./sqrt(phat.*(1-phat).*(1./n1 + 1./n2));
    end
end

pvals = 2*(1-normcdf(abs(zStats)));

% switch to 0/1/2/3 based on our criteria
zStats = abs(zStats);
for i = 1:numel(zStats)
     if zStats(i) > crit(3)
         zStats(i) = 3;
     elseif zStats(i) > crit(2)
         zStats(i) = 2;
     elseif zStats(i) >= crit(1)
         zStats(i) = 1;
     else
         zStats(i) = 0;
     end
end



end