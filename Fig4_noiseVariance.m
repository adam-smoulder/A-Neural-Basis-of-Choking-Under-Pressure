%%% This script looks at neural noise (i.e. trial-to-trial variability) as
%%% a function of reward. This script specifically looks at noise in the
%%% factor space.
%%%
%%% We'll look at variability in the following spaces:
%%% - full factor space
%%% - reward axis
%%% - target plane
%%% - target preparation axis
%%%
%%% Adam Smoudler, 12/15/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;

%% For each animal, get:
% - each noise variance listed in the header
% - error bars for them based on bootstrapping
method_R = 'PCA';
method_D = 'PCA';
orthogonalize = true;
nboots = 5000; % number of bootstraps for noise var error bars

% noiseVarLabels = {'Total','Reward axis','Target plane','TPA','Orthogonal'};
noiseVarLabels = {'Total','Reward axis','Target axes','Orthogonal'};

noiseVar_byRew_byCat_subjSplit = cell(nsubjects,1);
noiseVarErr_byRew_byCat_subjSplit = cell(nsubjects,1);
TPMean_byDir_byRew_subjSplit = cell(nsubjects,1);
TPCov_byDir_byRew_subjSplit = cell(nsubjects,1);
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
    
    % For variability, we don't want to include skip labels at all, as
    % outliers will dramatically skew calculations
    ntrials = sum(~skipLabels);
    directionLabels(skipLabels) = [];
    rewardLabels(skipLabels) = [];
    dayLabels(skipLabels) = [];
    delayLengths(skipLabels) = [];
    binnedData_pGC(skipLabels,:) = [];
    
    % Get target plane and projection
    [wD,muD,targPlaneProj,pveD] = getTargetPlane(binnedData_pGC,...
        rewardLabels,directionLabels,method_D);
    
    % Get noise covariance ellipses in the target plane
    for d = 1:ndirections
        for r = 1:nrewards
            curInds = directionLabels==directions(d) & rewardLabels==rewards(r);
            TPMean_byDir_byRew_subjSplit{f}(d,r,:) = mean(targPlaneProj(curInds,:));
            TPCov_byDir_byRew_subjSplit{f}(d,r,:,:) = cov(targPlaneProj(curInds,:));
        end; clear r
    end; clear d

    % Get reward axis *orthogonal to target plane*
    [wRnD,~,rewAxisProj,pveR] = getRewardAxis(binnedData_pGC*null(wD'),...
        rewardLabels,directionLabels,method_R);
    
    % And finally, get data orthogonal to the reward/target stuff
    orthData = binnedData_pGC*null(wD')*null(wRnD');
    
    % Get indices for each reward/direction condition
    [~,inds_byDir_byRew,n_byDir_byRew] = ...
        groupDataByLabel(directionLabels,[directionLabels rewardLabels]);
    
    % First, get the noise variance for each data subset
    varTotal_byRew = mean(cellfun(@(x) sum(var(binnedData_pGC(x,:))), inds_byDir_byRew));
    varTargPlane_byRew = mean(cellfun(@(x) sum(var(targPlaneProj(x,:))), inds_byDir_byRew));
    varRewAxis_byRew = mean(cellfun(@(x) var(rewAxisProj(x)), inds_byDir_byRew));
    varOrth_byRew = mean(cellfun(@(x) sum(var(orthData(x,:))), inds_byDir_byRew));
    
    % Now, do bootstrapping to get errorbars on this
    clear varTotal_byRew_byBoot varTargPlane_byRew_byBoot varTPR_byRew_byBoot varOffTPR_byRew_byBoot varRewAxis_byRew_byBoot
    for b = 1:nboots
        curBootInds = cellfun(@(x) x(randi(length(x),[length(x),1])),inds_byDir_byRew,'UniformOutput',false);
        varTotal_byRew_byBoot(:,b) = mean(cellfun(@(x) sum(var(binnedData_pGC(x,:))), curBootInds));
        varTargPlane_byRew_byBoot(:,b) = mean(cellfun(@(x) sum(var(targPlaneProj(x,:))), curBootInds));
        varRewAxis_byRew_byBoot(:,b) = mean(cellfun(@(x) var(rewAxisProj(x)), curBootInds));
        varOrth_byRew_byBoot(:,b) = mean(cellfun(@(x) sum(var(orthData(x,:))), curBootInds));
    end; clear b
    varErrTotal_byRew = std(varTotal_byRew_byBoot');
    varErrTargPlane_byRew = std(varTargPlane_byRew_byBoot');
    varErrRewAxis_byRew = std(varRewAxis_byRew_byBoot');
    varErrOrth_byRew = std(varOrth_byRew_byBoot');

    % Store values
    noiseVar_byRew_byCat_subjSplit{f} = [varTotal_byRew' varRewAxis_byRew' varTargPlane_byRew' varOrth_byRew'];
    % noiseVar_byRew_byCat_subjSplit{f} = [mean(varTotal_byRew_byBoot,2) mean(varRewAxis_byRew_byBoot,2) mean(varTargPlane_byRew_byBoot,2) mean(varOrth_byRew_byBoot,2)];
    noiseVarErr_byRew_byCat_subjSplit{f} = [varErrTotal_byRew' varErrRewAxis_byRew' varErrTargPlane_byRew' varErrOrth_byRew'];
end; clear f


%% Show the noise covariance ellipses in the target plane
figure;

rewColors = getDistinctColors('SELECT_ORDER',8);
dirColors = getDistinctColors('SELECT_ORDER',16);
ms = 6;
lw = 0.5;
axis_subjSplit = {[-90 120 -110 100], [-115 115 -115 115], [-210 250 -265 195]};
subplotOrder_subjSplit = {[1 2 7 8], [3 4 9 10], [5 6 11 12]};
for f = 1:nsubjects
    for r = 1:nrewards
        subplot(2,2*nsubjects,subplotOrder_subjSplit{f}(r)); hold on
        curMeans = squeeze(TPMean_byDir_byRew_subjSplit{f}(:,r,:));
        curCovs = squeeze(TPCov_byDir_byRew_subjSplit{f}(:,r,:,:));
        
        % Plot the ring
        ringPoints = [curMeans ; curMeans(1,:)];
        plot(ringPoints(:,1),ringPoints(:,2),'-','linewidth',lw,...
            'color',rewColors{rewards(r)})
        
        % Plot cov ellipses
        if strcmp(subjectNames{f},'Prez')
            directions = 2:2:8;
        else
            directions = 1:8;
        end
        ndirections = length(directions);
        for d = 1:ndirections
            plotCovEllipse(curMeans(d,:),squeeze(curCovs(d,:,:)),rewColors{r}); % default is 1 SD in 2D
        end; clear d
        for d = 1:ndirections
            plot(curMeans(d,1),curMeans(d,2),'o','markersize',ms,'linewidth',lw,...
                'markeredgecolor',dirColors{directions(d)},'markerfacecolor',dirColors{directions(d)});
        end; clear d
        det(cov(curMeans(:,1:2),1))
        %     axis equal
        %     axis(flat([axisMin(1:2) axisMax(1:2)]')')
        axis(axis_subjSplit{f})
%         axis equal
        set(gca,'fontsize',14);
        set(gca,'fontname','arial')
        if r ==1 && f == 1
            xlabel('Targ. Ax. 1')
            ylabel('Targ. Ax. 2')
        end
        xticks([])
        yticks([])
    end; clear r
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 1510 416])

% Save it
figname = ['Fig4_noiseVar_covEllipses_targMethod' method_D]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

% %% Make figures showing these
% 
% figure
% lw = 2;
% labelOrder = [3 2 4 1]; % TP, rew, orth, total
% for f = 1:nsubjects
%     for L = 1:length(labelOrder)
%       subplot(length(noiseVarLabels),nsubjects,f+(L-1)*nsubjects)
%       ind = labelOrder(L);
%       errorbar(1:nrewards,noiseVar_byRew_byCat_subjSplit{f}(:,ind),noiseVarErr_byRew_byCat_subjSplit{f}(:,ind),...
%           'linewidth',lw,'color',[0.8 0 0.8])
%       axis([0.5 4.5 -inf inf])
%       set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
%       xticks(rewards)
%       xticklabels({'S','M','L','J'})
%       if f == 1
%           ylabel([noiseVarLabels{ind}])
%       end
%       if L==1
%           title(['Monkey ' subjectNames{f}(1)])
%       end
%     end; clear L
% end; clear f
% pos = get(gcf,'position');
% set(gcf,'position',[pos(1:2)-[0 500] 813 867])
% 
% 
% % Save it
% figname = ['Fig4_noiseVar_rewMethod' method_R '_targMethod' method_D]
% saveFigAndSvg([figfolder 'recentFigs\'],figname);
% saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

%%


figure; set(gcf,'position',[-1588 70 1189 867])
lw = 2;
labelOrder = [1 3 2 4]; % TP, rew, orth, total
for f = 1:nsubjects
    for L = 1:length(labelOrder)
      subplot(nsubjects,length(noiseVarLabels),L+(f-1)*length(labelOrder))
      ind = labelOrder(L);
      errorbar(1:nrewards,noiseVar_byRew_byCat_subjSplit{f}(:,ind),noiseVarErr_byRew_byCat_subjSplit{f}(:,ind),...
          'linewidth',lw,'color',[0.8 0 0.8])
      axis([0.5 4.5 -inf inf])
      set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
      xticks(rewards)
      xticklabels({'S','M','L','J'})
      if L==1
          ylabel(['Monkey ' subjectNames{f}(1)])
      end
      if f==1
          title([noiseVarLabels{ind}])
      end
    end; clear L
end; clear f

% Save it
figname = ['Fig4_noiseVar_rewMethod' method_R '_targMethod' method_D]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);