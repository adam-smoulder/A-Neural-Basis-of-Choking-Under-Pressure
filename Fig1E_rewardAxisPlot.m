%%% This script makes plots showing the reward axis found using PCA on
%%% reward means, agnostic to target condition.
%%%
%%% Adam Smoulder, 7/14/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;


%% For each animal, get the reward axis and projections along it, split by 
%  reward
rewMethod = 'PCA';
subjectNames = cell(nsubjects,1);
varExp_subjSplit = cell(nsubjects,1);
varExpNull_subjSplit = cell(nsubjects,1);
rewProj_subjSplit = cell(nsubjects,1);
rewardLabels_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    nfactors = size(binnedData_pGC,2);
    
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
    
    n_subjSplit{f} = nByLabel(rewardLabels);
    
    % Get reward axis
    [wR,muR,~,varExp_subjSplit{f}] = getRewardAxis(binnedData_pGC(~skipLabels,:),...
        rewardLabels(~skipLabels,:),directionLabels(~skipLabels,:),rewMethod);

    % Project all data down to the reward axis store
    rewProjData = (binnedData_pGC-muR)*wR;
    rewProj_subjSplit{f} = rewProjData(~skipLabels);
    rewardLabels_subjSplit{f} = rewardLabels(~skipLabels);
    
    % Get null boundary for var exp; we'll randomly pick 4 points from the
    % covariance of the data and do PCA on those like we did with reward
    % means, then make a null dist of what we'd expect varExp for PC1 to be
    nreps = 10000;
    sigma = cov(binnedData_pGC(~skipLabels,:));
    varExpNull_byRep = nan(nreps,3);
    for n = 1:nreps
        curNullData = mvnrnd(zeros(nfactors,1),sigma,4);
        [~,~,~,~,varExpNull_byRep(n,:)] = pca(curNullData-mean(curNullData));
    end
    varExpNull_subjSplit{f} = varExpNull_byRep;
    
end; clear f
disp('Done processing data')

% We'll call "significantly 1D" as PC1 explaining more variance than 95% of
% the null distribution. The p-value then is 1-(what percentile of the null
% the data is at)
oneDpValues = 1-cellfun(@(x,y) sum(x(1) > y(:,1))/length(y), varExp_subjSplit, varExpNull_subjSplit)


%% Plot medians and individual data points
if strcmp(rewMethod,'PCA') % use the pre-set vals
    minVal_bySubj = [-120 -120 -275];
    maxVal_bySubj = [100 200 225];
else % calculate
    minVal_bySubj = nan(nsubjects,1);
    maxVal_bySubj = nan(nsubjects,1);
    for f = 1:nsubjects
        minVal_bySubj(f) = 1.1*prctile(rewProj_subjSplit{f},0.5);
        maxVal_bySubj(f) = 1.1*prctile(rewProj_subjSplit{f},99.5);
    end; clear f
end


rewNames = {'Small','Medium','Large','Jackpot'};
rewColors = getDistinctColors('SELECT_ORDER',8);
colorRewNames = rewNames;
for r = 1:nrewards
    colorRewNames{r} = ['\color{red} ' rewNames{r}];
    colorRewNames{r} = sprintf('\\color[rgb]{%f,%f,%f}%s', 0.9*rewColors{r}, rewNames{r});
end; clear r

figure
lw = 2;
ms1 = 5;
ms2 = 8;
alphas_subjSplit = {...
    0.25*[0.3 0.5 0.3 0.5],...
    0.3*[0.3 0.5 0.3 0.5],...
    0.15*[0.3 0.5 0.3 0.5]
    };
for f = 1:nsubjects
    subplot(1,nsubjects,f); hold on
    rewProj_byRew = groupDataByLabel(rewProj_subjSplit{f},rewardLabels_subjSplit{f});
    rewProj_byRew_median = cellfun(@(x) median(x), rewProj_byRew);
    rewProj_byRew_semed = cellfun(@(x) nansemed(x,1), rewProj_byRew);
    
    % Set up axes
    minVal = minVal_bySubj(f);
    maxVal = maxVal_bySubj(f);
    
    % Plot datapts
    for r = 1:nrewards
        y = rewProj_byRew{r};
        x = r*ones(size(y))+0.15*randn(size(y));
        scatter(x,y,ms1,rewColors{r},'filled','markerfacealpha',alphas_subjSplit{f}(r),'markeredgealpha',alphas_subjSplit{f}(r))
    end; clear r
    
    % Plot the medians
    plot(1:nrewards,rewProj_byRew_median,'k-','linewidth',lw)
    for r = 1:nrewards
        plot(r,rewProj_byRew_median(r),'o','linewidth',lw,'markersize',ms2,...
            'color',0.7*rewColors{r},'markerfacecolor',rewColors{r})
    end; clear r
    
    % Labels and such
    if f == 1
        ylabel('Reward Axis (Hz)')
    end
%     if f == 2
%         xlabel('Reward')
%     end
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    title(['Monkey ' subjectNames{f}(1)])
    xticks(1:4)
    xticklabels(colorRewNames)
        xtickangle(45)

    axis([0.25 4.75 minVal maxVal])
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 724 432])

% Save it
figname = ['Fig1D_rewardAxisPlot_' rewMethod ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);




% 
% %% Make violin plots of the data
% 
% figure
% rewColors = getDistinctColors('SELECT_ORDER',8);
% lw = 2;
% for f = 1:nsubjects
%     subplot(1,nsubjects,f); hold on
%     rewProj_byRew = groupDataByLabel(rewProj_subjSplit{f},rewardLabels_subjSplit{f});
%     rewProj_byRew_mean = cellfun(@(x) mean(x), rewProj_byRew);
%     rewProj_byRew_sem = cellfun(@(x) nansem(x,1), rewProj_byRew);
%     
%     % Set up axes and histograms
% %     minVal = prctile(rewProj_subjSplit{f},0.1);
% %     maxVal = 1.2*prctile(rewProj_subjSplit{f},99.99);
%     switch subjectNames{f}
%         case 'Earl'
%             minVal = -125; maxVal = 125;
%         case 'Prez'
%             minVal = -155; maxVal = 255;
%         case 'Rocky'
%             minVal = -305; maxVal = 255;
%     end
%     nbins = 50;
%     histBinEdges = linspace(minVal,maxVal,nbins+1);
%     histBinCenters = 1/2*(histBinEdges(1:end-1)+histBinEdges(2:end));
%     rewProjPmfs_byRew = nan(nrewards,nbins);
%     for r = 1:nrewards
%         rewProjPmfs_byRew(r,:) = histcounts(rewProj_byRew{r},histBinEdges,'normalization','probability');
%     end; clear r
%     
%     % Smooth the pmf lightly!!! Not too hard. Just enough to make violin-like
%     kernelStd = 1.25;
%     kernel = normpdf(-6*kernelStd:6*kernelStd,0,kernelStd);
%     rewProjPmfs_byRew_smooth = convByDim(rewProjPmfs_byRew,kernel,2,'same');
%     rewProjPmfs_byRew_smooth((rewProjPmfs_byRew_smooth < 3E-3)) = nan;
% 
%     % Plot the violins
%     for r = 1:nrewards
%         curPmf = rewProjPmfs_byRew_smooth(r,:);
%         curPmf_norm = curPmf./max(rewProjPmfs_byRew_smooth(:))*0.4; % scale all equally
%         x1 = r-curPmf_norm;
%         x2 = r+curPmf_norm;
%         y = histBinCenters;
%         plot(x1,y,'-','linewidth',lw,'color',rewColors{r});
%         plot(x2,y,'-','linewidth',lw,'color',rewColors{r});
%         plot(r*ones(100,1),linspace(minVal,maxVal,100),'-','color',rewColors{r},'linewidth',0.25)
%     end; clear r
%     plot(1:nrewards,rewProj_byRew_mean,'k-','linewidth',1)
%     for r = 1:nrewards
% %         errorbar(rewards(r),rewProj_byRew_mean(r),rewProj_byRew_sem(r),'o',...
% %             'color',rewColors{r},'markerfacecolor',rewColors{r},'markersize',2,'linewidth',1.5,'capsize',8)
%         plot(r,rewProj_byRew_mean(r),'o','color',rewColors{r},'markerfacecolor',rewColors{r},'markersize',8)
%     end; clear r
%     
%     % Labels and such
%     if f == 1
%         ylabel('Reward Axis (Hz)')
%     end
%     if f == 2
%         xlabel('Reward')
%     end
%     set(gca,'fontsize',12,'fontname','arial','tickdir','out');
%     title(['Monkey ' subjectNames{f}(1)])
%     xticks(1:4)
%     xticklabels({'S','M','L','J'})
% %     yticks([minVal+5 maxVal-5])
%     yticks([-50 0 50])
%     if f~=1
%         yticklabels({[],'0',[]})
%     end
%     axis([0.25 4.75 minVal maxVal])
% end; clear f
% set(gcf,'position',[504 511 710 437])
% 
% % Save it
% if strcmp(rewMethod,'PCA')
%     figname = ['Fig1D_rewardAxisViolins_' rewMethod ]
% else
%     figname = ['FigS3_rewardAxisViolins_' rewMethod ]
% end
% saveFigAndSvg([figfolder 'recentFigs\'],figname);
% saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
