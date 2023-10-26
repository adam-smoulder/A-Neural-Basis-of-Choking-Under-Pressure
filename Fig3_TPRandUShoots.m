%%% This script looks at the target plane radius for successes versus
%%% overshoots vs. undershoots.
%%%
%%% Adam Smoulder, 8/29/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;


%% Get target plane radius and correlate with behaviors
method = 'PCA_dirZ'; % 'PCA'
ndimout = 2;
orthogonalize = true;
subjectNames = cell(nsubjects,1);
TPRProj_subjSplit = cell(nsubjects,1);
directionLabels_subjSplit = cell(nsubjects,1);
rewardLabels_subjSplit = cell(nsubjects,1);
statusLabels_subjSplit = cell(nsubjects,1);
dayLabels_subjSplit = cell(nsubjects,1);
prewLabels_subjSplit = cell(nsubjects,1);
TPRProj_byStatus_byRew_mean_subjSplit = cell(nsubjects,1);
TPRProjRewZ_byStatus_subjSplit = cell(nsubjects,1);
nrewards = 4;
pvals_oshoot_bySubj_byRew = nan(nsubjects,nrewards);
pvals_ushoot_bySubj_byRew = nan(nsubjects,nrewards);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    [ntrials,nfactors] = size(binnedData_pGC);
    
    % Remove outliers or short delay trials as needed
    if strcmp(taskInfo.subjectName,'Earl') % Earl has some 
        skipLabels = sum(isoutlier(binnedData_pGC,'thresholdfactor',4),2)>0; % very outlier-y, but present enough to influence PCA
    elseif strcmp(taskInfo.subjectName,'Rocky')
        skipLabels = delayLengths < 317; % excluding trials w/ delay shorter than dynamics settling
        % skipLabels = delayLengths <= 449; % static value to match monkey E's
    else % Prez's neural activity asymptotes in time, so no bad inds
        skipLabels = false(length(binnedData_pGC),1);
        % skipLabels = delayLengths <= 449; % static value to match monkey E's
    end
%     skipLabels = false(length(binnedData_pGC),1); % uncomment to include all trials for all animals
    disp(['n trials to skip = ' num2str(sum(skipLabels))])
    
    % Get target plane radius and projection
    if ~contains(method,'full') % we do not need to get target plane if using full data
        [wD,muD] = getTargetPlane(binnedData_pGC(~skipLabels,:),...
            rewardLabels(~skipLabels,:),directionLabels(~skipLabels,:),method);
        targPlaneProj = (binnedData_pGC-muD)*wD;
    else % we'll use the full data for the target plane
        targPlaneProj = binnedData_pGC;
    end
    
    % Now get target plane radius by projecting trials onto the center of
    % the vector connecting the center of their condition's reward ring to
    % their direction average
    TPRProjData = nan(ntrials,1);
    for d = 1:ndirections
        curInds = directionLabels==directions(d);
        targPlaneProj_byDir_byRew_mean = cellfun(@(x) mean(x),...
            groupDataByLabel(targPlaneProj(~skipLabels,:),...
            [directionLabels(~skipLabels) rewardLabels(~skipLabels)]),'uniformoutput',false);
        center_byRew = squeeze(mean(...
            reshape([targPlaneProj_byDir_byRew_mean{:}],...
            [size(targPlaneProj,2) ndirections nrewards]),2))';
        for r = 1:nrewards
            curInds = directionLabels==directions(d) & rewardLabels==rewards(r);
            endpt = targPlaneProj_byDir_byRew_mean{d,r};
            unitVec = endpt-center_byRew(r,:);
            unitVec = unitVec/norm(unitVec);
            TPRProjData(curInds) = (targPlaneProj(curInds,:)-center_byRew(r,:))*unitVec';
        end; clear r
    end; clear d
    
    % Z-score by direction if desired
    if contains(method,'dirZ')
        TPRProj_byDir = groupDataByLabel(TPRProjData(~skipLabels),[directionLabels(~skipLabels)]);
        ndirections = size(TPRProj_byDir,1); directions = unique(directionLabels);
        mu_byDir = nan(ndirections,1);
        sigma_byDir = nan(ndirections,1);
        for d = 1:ndirections
            curInds = directionLabels==directions(d);
            curZInds = directionLabels==directions(d) & trialStatusLabels==1;
            mu_byDir(d) = mean(TPRProjData(curZInds));
            sigma_byDir(d) = std(TPRProjData(curZInds));
            TPRProjData(curInds) = (TPRProjData(curInds)-mu_byDir(d))/sigma_byDir(d);
        end; clear d
        
%         % If we want it to not be centered at zero bc bars look worse, use
%         % the average value scaled by the avg std
%         TPRProjData = TPRProjData+mean(mu_byDir./sigma_byDir);
    end
    TPRProj_subjSplit{f} = TPRProjData(~skipLabels);
    directionLabels_subjSplit{f} = directionLabels(~skipLabels);
    rewardLabels_subjSplit{f} = rewardLabels(~skipLabels);
    statusLabels = [... % 1 = succ, 2 = oshoot, 3 = ushoot
        ismember(trialStatusLabels,1)...
        ismember(trialStatusLabels,[-22 -34])...
        ismember(trialStatusLabels,-23)]*[1 2 3]';
    statusLabels_subjSplit{f} = statusLabels(~skipLabels);
    dayLabels_subjSplit{f} = dayLabels(~skipLabels);
    prewLabels_subjSplit{f} = prevRewardLabels(~skipLabels);
    
    % Split by success and failure; save avgs and zs by reward
    TPRProj_byStatus_byRew = groupDataByLabel(TPRProjData(~skipLabels),[statusLabels(~skipLabels), rewardLabels(~skipLabels)]);
    n_byStatus_byRew_subjSplit{f} = cellfun(@(x) length(x), TPRProj_byStatus_byRew);
    TPRProj_byStatus_byRew_mean = cellfun(@(x) mean(x), TPRProj_byStatus_byRew);
    TPRProj_byStatus_byRew_sem = cellfun(@(x) nansem(x,1), TPRProj_byStatus_byRew);
    TPRProj_byStatus_byRew_mean_subjSplit{f} = TPRProj_byStatus_byRew_mean;
    TPRProj_byStatus_byRew_sem_subjSplit{f} = TPRProj_byStatus_byRew_sem;
    
    
    % Z-score within reward
    TPRProjRewZ = nan(size(TPRProjData));
    for r = 1:nrewards % z-score by reward success
        curInds = rewardLabels==rewards(r) & ~skipLabels;
        curProjData = TPRProjData(curInds);
        curZData = TPRProjData(curInds & trialStatusLabels==1);
        TPRProjRewZ(curInds) = (curProjData-mean(curZData))./std(curZData);
    end; clear r
    TPRProjRewZ_byStatus = groupDataByLabel(TPRProjRewZ(~skipLabels),statusLabels(~skipLabels));
    TPRProjRewZ_byStatus_subjSplit{f} = TPRProjRewZ_byStatus;
    
%     % There's some weird issues with the ttest? Let's just use
%     % bootstrapping
    [~,pvals_oshoot_bySubj_byRew_tt(f,:)] = cellfun(@(x,y) ttest2(x,y),TPRProj_byStatus_byRew(1,:),TPRProj_byStatus_byRew(2,:));
    [~,pvals_ushoot_bySubj_byRew_tt(f,:)] = cellfun(@(x,y) ttest2(x,y),TPRProj_byStatus_byRew(1,:),TPRProj_byStatus_byRew(3,:));
    
    % We'll do a bootstrap test to determine a sig. difference in means;
    % see that it corroborates ttest
    nboots = 1E5;
    for r = 1:nrewards
        succ_oshoot_bootMeans = nan(nboots,1);
        succ_ushoot_bootMeans = nan(nboots,1);
        x = TPRProj_byStatus_byRew{1,r};
        y = TPRProj_byStatus_byRew{2,r};
        z = TPRProj_byStatus_byRew{3,r};
        for b = 1:nboots
            xInds = randi(length(x),[length(x) 1]);
            succ_oshoot_bootMeans(b) = mean(x(xInds))-mean(y(randi(length(y),[length(y) 1])));
            succ_ushoot_bootMeans(b) = mean(x(xInds))-mean(z(randi(length(z),[length(z) 1])));
        end; clear b
        pvals_oshoot_bySubj_byRew(f,r) = min(sum(succ_oshoot_bootMeans > 0),sum(succ_oshoot_bootMeans < 0))/nboots*2; % we check whichever is lower and mult x 2 bc it's 2-sided
        pvals_ushoot_bySubj_byRew(f,r) = min(sum(succ_ushoot_bootMeans > 0),sum(succ_ushoot_bootMeans < 0))/nboots*2; % we check whichever is lower and mult x 2 bc it's 2-sided
    end; clear r
    

end; clear f
disp('Done processing data')


%% Make the slope plot

figure
subjShades = [1 0.75 0.5];
subjShapes = 'os^';
rewColors = getDistinctColors('SELECT_ORDER',8);
lw = 2.5;
ms = 5;
minVal = min(cellArrayToVector(TPRProj_byStatus_byRew_mean_subjSplit));
maxVal = max(cellArrayToVector(TPRProj_byStatus_byRew_mean_subjSplit));
maxVal = maxVal+0.25*(maxVal-minVal);
minVal = minVal-0.15*(maxVal-minVal);

subplot(1,2,2); hold on % show success vs. overshoot

% Plot a fake value out of range for the legend's sake
for f = 1:nsubjects
    plot(1,maxVal*10,subjShapes(f),'markeredgecolor','k',...
        'linewidth',lw,'markersize',ms)
end; clear f

% Now plot the real stuff (overshoots)
for f = 1:nsubjects
    xmod = -0.15+f*0.075;
    for r = 1:nrewards
        curVals = TPRProj_byStatus_byRew_mean_subjSplit{f}(:,r);
%         curVals = curVals-curVals(1); % set success to 0
        plot((1:2)+xmod,curVals(1:2),'k-','linewidth',0.5+(pvals_oshoot_bySubj_byRew(f,r)<0.05)*2,'color',[0.7 0.7 0.7])
        plot(1+xmod,curVals(1),subjShapes(f),'markeredgecolor',rewColors{r},...
            'linewidth',lw,'markersize',ms)
        plot(2+xmod,curVals(2),subjShapes(f),'markeredgecolor',rewColors{r},...
            'linewidth',lw,'markersize',ms)
    end; clear r
end; clear f
% legend(cellfun(@(x) x(1),subjectNames,'uniformoutput',false),'location','NE')
axis([0.5 2.5 minVal maxVal])
xticks(1:2); xticklabels({'Success','Overshoot'}); xtickangle(22.5)
yticks([-1.5:0.5:1.5]); yticklabels([])
set(gca,'fontsize',12,'fontname','arial','tickdir','out');
legend(cellfun(@(x) x(1),subjectNames,'uniformoutput',false),'location','SE')



subplot(1,2,1); hold on
% Plot a fake value out of range for the legend's sake
for f = 1:nsubjects
    plot(1,maxVal*10,subjShapes(f),'markeredgecolor','k',...
        'linewidth',lw,'markersize',ms)
end; clear f

% Now plot the real stuff
for f = 1:nsubjects
    xmod = -0.15+f*0.075;
    for r = 1:nrewards
        curVals = TPRProj_byStatus_byRew_mean_subjSplit{f}(:,r);
%         curVals = curVals-curVals(1); % set success to 0
        plot((1:2)+xmod,curVals([1 3]),'k-','linewidth',0.5+(pvals_ushoot_bySubj_byRew(f,r)<0.05)*2,'color',[0.7 0.7 0.7])
        plot(1+xmod,curVals(1),subjShapes(f),'markeredgecolor',rewColors{r},...
            'linewidth',lw,'markersize',ms)
        plot(2+xmod,curVals(3),subjShapes(f),'markeredgecolor',rewColors{r},...
            'linewidth',lw,'markersize',ms)
    end; clear r
end; clear f
axis([0.5 2.5 minVal maxVal])
xticks(1:2); xticklabels({'Success','Undershoot'}); xtickangle(22.5)
yticks([-1.5:0.5:1.5])
set(gca,'fontsize',12,'fontname','arial','tickdir','out');
ylabel('Target Preparation Axis')

set(gcf,'position',[560 528 450 420])

% Save
figname = ['Fig3_TPRandUShoots_rewSFCombinedSubj_' method]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


% Show t-test difference between two
data = [TPRProj_byStatus_byRew_mean_subjSplit{:}]'; % [datapts x S/F]
[~,p_oshoot,ci_oshoot,stats_oshoot] = ttest(data(:,1),data(:,2));
[~,p_ushoot,ci_ushoot,stats_ushoot] = ttest(data(:,1),data(:,3));
disp(['p (succ/oshoot) = ' num2str(p_oshoot)])
disp(['p (succ/ushoot) = ' num2str(p_ushoot)])

%% For each subject, calculate the error bars for the difference b/w success
% and failure, then plot them
figure
meanVals = nan(nsubjects,2);
semVals = nan(nsubjects,2);
pVals_oshoot = nan(nsubjects,1);
pVals_ushoot = nan(nsubjects,1);
for f = 1:nsubjects
    curVals_byStatus = TPRProjRewZ_byStatus_subjSplit{f};
    meanVals(f,:) = cellfun(@(x) mean(x)-mean(curVals_byStatus{1}),curVals_byStatus(2:3));
    semVals(f,1) = sqrt(sum(cellfun(@(x) var(x)/length(x),curVals_byStatus(1:2))));
    semVals(f,2) = sqrt(sum(cellfun(@(x) var(x)/length(x),curVals_byStatus([1 3]))));
    [~,pVals_oshoot(f)] = ttest2(curVals_byStatus{1},curVals_byStatus{2},'vartype','unequal');
    [~,pVals_ushoot(f)] = ttest2(curVals_byStatus{1},curVals_byStatus{3},'vartype','unequal');
end; clear f
disp(['Pvals oshoot = ' num2str(pVals_oshoot')])
disp(['Pvals ushoot = ' num2str(pVals_ushoot')])

figure; hold on
lw = 1;
statusColors = {[0 1 0],[0.3 0.3 0.3],[0.7 0.7 0.7]};
barWidth = 0.38;
bar((1:nsubjects)+barWidth/2,meanVals(:,1),'facecolor',statusColors{2},'facealpha',0.5,'barwidth',barWidth)
errorbar((1:nsubjects)+barWidth/2,meanVals(:,1),semVals(:,1),'k.','linewidth',lw)
bar((1:nsubjects)-barWidth/2,meanVals(:,2),'facecolor',statusColors{3},'facealpha',0.5,'barwidth',barWidth)
errorbar((1:nsubjects)-barWidth/2,meanVals(:,2),semVals(:,2),'k.','linewidth',lw)

% Label stuff
axis([0.5 3.5 min(meanVals(:)-semVals(:))*1.1 0.14])
% yticks([-0.6:0.05:0.1])
set(gca,'fontsize',12,'fontname','arial','tickdir','out');
ylabel('\Delta TPA_{succ,fail} (z-score)')
xticks(1:nsubjects)
xticklabels(cellfun(@(x) x(1),subjectNames,'uniformoutput',false))
set(gca,'xaxislocation','top')
set(gcf,'position',[560 528 453 420])


% Save 
figname = ['Fig3_TPRAndUShoots_SFRewZBar_' method]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
