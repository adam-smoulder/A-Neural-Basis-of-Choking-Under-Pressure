function [decodeAcc_mean,decodeAcc_lowerErrorBar,decodeAcc_upperErrorBar]...
    = decodeAccPlot(decodeAccs,Nr)
% This function plots decode accuracy with a 95% CI and evalutates
% significant difference between different reward size decoders. The input
% is the bootstrapped distributions for each decoder that have been
% bias-corrected (means should be taken from subsampled decode accuracy to
% ensure more testing data isn't biasing decoders that have more samples).
% Significance is evaluated if 95% CIs of one distribution do not overlap
% the mean of another and indicate a difference in means between groups
% (rewards). This is equivalent to a t-test of sorts between the
% bootstrapped distributions. 
%
% Inputs:
% - decodeAccs: [N x K] matrix with N bootstrapped samples for each of K
% groups.
% - Nr:  minimum number of trials for a given cond size across directions
% (so how many pts go into the bootstraps), should be K x 1
%
% Outputs:
% - decodeAcc_mean: [K x 1] vector of mean decode accuracies
% - decodeAcc_lowerErrorBar: [K x 1] distance from the mean of lower 95% CI
% Bar
% - decodeAcc_upperErrorBar: [K x 1] same^ but for upper 95% CI bar
%
% Adam Smoulder, 8/8/19 (edit 1/20/20 for outputs)


alphaDiv2 = 0.025;
BC = 1; % Bonferonni correction

tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution


[nboots, nrewards] = size(decodeAccs);

% convert to percentages for ease of visualization
decodeAccs = 100*decodeAccs;

% evaluate nonoverlapping confidence intervals within each latent
decodeAcc_mean = squeeze(mean(decodeAccs));
acc2p5per = squeeze(prctile(decodeAccs,alphaDiv2/BC*100));
acc97p5per = squeeze(prctile(decodeAccs,(1-alphaDiv2/BC)*100));
allAccPs = [];
allAccRels = [];

for j = 1:(nrewards-1)
    for k = (j+1):nrewards
        accDiff = squeeze(decodeAccs(:,j)-decodeAccs(:,k));
        accDenom = sqrt(squeeze(var(decodeAccs(:,j))+var(decodeAccs(:,k))));
        accP = tdist2T(mean(accDiff)/accDenom,Nr(j)+Nr(k)-1);
        %             if (sum(snrDiff > 0)/length(snrDiff) > 0.95) || (sum(snrDiff < 0)/length(snrDiff) > 0.95)
        if accP <= alphaDiv2/BC || accP >= (1-alphaDiv2/BC)
            allAccRels = [allAccRels, {[j,k]}];
            allAccPs(end+1) = accP;
        end
    end; clear k
end; clear j


% plot it!
colors = getDistinctColors;
figure
hold on
plot(1:4,decodeAcc_mean(:),'k-','linewidth',2)
decodeAcc_lowerErrorBar = acc2p5per-decodeAcc_mean;
decodeAcc_upperErrorBar = acc97p5per-decodeAcc_mean;
for r = 1:4
    errorbar(r,decodeAcc_mean(r),decodeAcc_lowerErrorBar(r),decodeAcc_upperErrorBar(r),'ko','linewidth',3,'color',colors{r},'markerFaceColor',colors{r},'markeredgecolor',colors{r})
end
if ~isempty(allAccRels)
    sigstar(allAccRels, ternop(allAccPs>0.5,1-allAccPs,allAccPs));
end
xticks(1:4)
xticklabels({'S','M','L','J'});
ylabel('Decode Accuracy (%)')
title('Naive Bayes Decode Performance')
set(gca,'fontsize',18)
axis([0.8 4.2 -inf inf])

end

