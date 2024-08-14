%%% This script gets the target plane radius as a function of reward for
%%% each animal.
%%%
%%% Adam Smoulder, 7/19/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;

%% For each animal, get the target plane
ndimout = 2;
method = 'PCA_dirZ'; %'full_dirZ'
orthogonalize = true
subjectNames = cell(nsubjects,1);
R2_subjSplit = cell(nsubjects,1);
TPRProj_subjSplit = cell(nsubjects,1);
directionLabels_subjSplit = cell(nsubjects,1);
rewardLabels_subjSplit = cell(nsubjects,1);
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
        targPlaneProj = (binnedData_pGC(~skipLabels,:)-muD)*wD;
        disp(['Angle between target plane dims (deg): ' num2str(acosd(wD(:,1)'*wD(:,2)/(norm(wD(:,1))*norm(wD(:,2)))))])
    else % we'll use the full data for the target plane
        targPlaneProj = binnedData_pGC(~skipLabels,:);
    end
    directionLabels_subjSplit{f} = directionLabels(~skipLabels);
    rewardLabels_subjSplit{f} = rewardLabels(~skipLabels);
        
    % Calculate target plane radius
    TPRProj_subjSplit{f} = nan(sum(~skipLabels),1);
    targPlaneProj_byDir_byRew_mean = cellfun(@(x) mean(x),...
        groupDataByLabel(targPlaneProj,[directionLabels(~skipLabels) rewardLabels(~skipLabels)]),'uniformoutput',false);
    center_byRew = squeeze(mean(...
        reshape([targPlaneProj_byDir_byRew_mean{:}],[size(targPlaneProj,2) ndirections nrewards]),2))';
    for d = 1:ndirections
        for r = 1:nrewards
            curInds = directionLabels(~skipLabels)==directions(d) & rewardLabels(~skipLabels)==rewards(r);
            endpt = targPlaneProj_byDir_byRew_mean{d,r};
            unitVec = endpt-center_byRew(r,:);
            unitVec = unitVec/norm(unitVec);
            TPRProj_subjSplit{f}(curInds) = (targPlaneProj(curInds,:)-center_byRew(r,:))*unitVec';
        end; clear r
    end; clear d
end; clear f
disp('Done processing data')


% Z-score within each direction
for f = 1:nsubjects
    if contains(method,'dirZ')
        TPRProj_byDir = groupDataByLabel(TPRProj_subjSplit{f},[directionLabels_subjSplit{f}]);
        ndirections = size(TPRProj_byDir,1); directions = unique(directionLabels_subjSplit{f});
        for d = 1:ndirections
            curInds = directionLabels_subjSplit{f}==directions(d);
%             TPRProj_subjSplit{f}(curInds) = (TPRProj_subjSplit{f}(curInds)-mean(TPRProj_byDir{d}))/std(TPRProj_byDir{d});
            TPRProj_subjSplit{f}(curInds) = (TPRProj_subjSplit{f}(curInds)-median(TPRProj_byDir{d}))/(0.741*iqr(TPRProj_byDir{d})); % robust z-score
        end; clear d
    end
end; clear f



%% Plot datapts and medians
figure
rewColors = getDistinctColors('SELECT_ORDER',8);
lw = 2;
ms1 = 5;
ms2 = 8;
alphas_subjSplit = {...
    0.25*[0.3 0.5 0.3 0.5],...
    0.3*[0.3 0.5 0.3 0.5],...
    0.15*[0.3 0.5 0.3 0.5]
    };
pvalMat_subjSplit = cell(nsubjects,1);
TPRProj_byRew_mean_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    subplot(1,nsubjects,f); hold on
    TPRProj_byRew = groupDataByLabel(TPRProj_subjSplit{f},rewardLabels_subjSplit{f});
%     TPRProj_byRew_val = cellfun(@(x) median(x), TPRProj_byRew);
%     TPRProj_byRew_err = cellfun(@(x) nansemed(x,1), TPRProj_byRew);
    TPRProj_byRew_val = cellfun(@(x) mean(x), TPRProj_byRew);
    TPRProj_byRew_err = cellfun(@(x) nansem(x,1), TPRProj_byRew);

    % Plot datapts
    for r = 1:nrewards
        y = TPRProj_byRew{r}; % success
        x = r*ones(size(y))+0.15*randn(size(y));
        scatter(x,y,ms1,rewColors{r},'filled','markerfacealpha',alphas_subjSplit{f}(r),'markeredgealpha',alphas_subjSplit{f}(r))
    end; clear r
    
    % Plot the bars
    plot(1:nrewards,TPRProj_byRew_val,'k-','linewidth',lw)
    for r = 1:nrewards
        plot(r,TPRProj_byRew_val(r),'o','linewidth',lw,'markersize',ms2,...
            'color',0.7*rewColors{r},'markerfacecolor',rewColors{r})
    end; clear r
    
    
    % Labels and such
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    title(['Monkey ' subjectNames{f}(1)])
    xticks(1:4)
    xticklabels({'S','M','L','J'})
    if contains(method,'dirZ')
        minVal = -3;
        maxVal = 3;
        yticks([-2.5 0 2.5])
        yaxistitle = 'Target Preparation Axis (a.u.)';
    else
        minVal = 1.1*prctile(TPRProj_subjSplit{f},0.5);
        maxVal = 1.1*prctile(TPRProj_subjSplit{f},99.5);
        yaxistitle = 'Target Projection Axis (sp/s)';
    end
    if f == 1
        ylabel(yaxistitle)
    end
    if f == 2
        xlabel('Reward')
    end
    axis([0.25 4.75 minVal maxVal])
    
    % Get p-values with Welch's t-test
    pvalMat_subjSplit{f} = nan(nrewards,nrewards);
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            [~,pvalMat_subjSplit{f}(r1,r2)] = ttest2(TPRProj_byRew{r1},TPRProj_byRew{r2});
        end; clear r2
    end; clear r1
    TPRProj_byRew_mean_subjSplit{f} = TPRProj_byRew_val;
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 710 437])

% Save it!
figname = ['Fig2C_targetPlaneRadiusByReward_datapts_' method '_orth' num2str(orthogonalize)]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

% Display p-values
for f = 1:nsubjects
    disp(subjectNames{f})
    disp([num2str(diff(TPRProj_byRew_mean_subjSplit{f}([1 3]))) ', p_{S->L} = ' num2str(pvalMat_subjSplit{f}(1,3))])
    disp([num2str(diff(TPRProj_byRew_mean_subjSplit{f}([3 4]))) ', p_{L->J} = ' num2str(pvalMat_subjSplit{f}(3,4))])
end; clear f

