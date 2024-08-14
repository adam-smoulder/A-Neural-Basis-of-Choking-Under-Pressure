%%% This script looks at the relationship between target plane radius and
%%% behavior. Also the "off-TPA" and behavior.
%%%
%%% Adam Smoulder, 4/17/24


addpath(genpath('C:\Users\Smoulgari\Documents\MATLAB\'))
figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
filename_bySubject = {...
    'NeuralDataForFigs_Earl_20240522_125842'; % Update inclues angular error
    'NeuralDataForFigs_Prez_20240522_125959';
    'NeuralDataForFigs_Rocky_20240522_130143';
    };
nsubjects = length(filename_bySubject);
dateString = grabDateTimeString;
rng(3195)

%% Get target plane radius and correlate with behaviors
method = 'PCA_dirZ'; %'full_dirZ'; % 'PCA'
ndimout = 2;
orthogonalize = true;
subjectNames = cell(nsubjects,1);
TPRProj_subjSplit = cell(nsubjects,1);
directionLabels_subjSplit = cell(nsubjects,1);
rewardLabels_subjSplit = cell(nsubjects,1);
successLabels_subjSplit = cell(nsubjects,1);
dayLabels_subjSplit = cell(nsubjects,1);
prewLabels_subjSplit = cell(nsubjects,1);
corr_RT_TPR_byDR_subjSplit = cell(nsubjects,1);
corr_PS_TPR_byDR_subjSplit = cell(nsubjects,1);
corr_HT_TPR_byDR_subjSplit = cell(nsubjects,1);
corr_BEP_TPR_byDR_subjSplit = cell(nsubjects,1);
% corr_AE_TPR_byDR_subjSplit = cell(nsubjects,1);
corr_absoffPSL_TPR_byDR_subjSplit = cell(nsubjects,1); % off-target axis peak speed location (e.g., off axis deviation)
corr_PSgRT_TPR_byDR_subjSplit = cell(nsubjects,1); % peak speed given reaction time (partial correlation)
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load([filename_bySubject{f}])
    subjectNames{f} = taskInfo.subjectName;
    [ntrials,nfactors] = size(binnedData_pGC);
    
    % Remove outlier trials before fitting a PCA model; Earl is the only
    % animal where this matters, as his F1 shows some crazy outliers
    if strcmp(taskInfo.subjectName,'Earl')
        skipLabels = sum(isoutlier(binnedData_pGC,'thresholdfactor',4),2)>0; % VERY outlier-y, but present enough to influence...
    elseif strcmp(taskInfo.subjectName,'Rocky')
        skipLabels = delayLengths < 317;
    else % Prez's neural activity asymptotes in time, so no bad inds
        skipLabels = false(length(binnedData_pGC),1);
    end
    disp(['n outliers/to skip = ' num2str(sum(skipLabels))])
    
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
    oTPRProjData = nan(ntrials,1); % off axis TPA
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
            rotMat = [cosd(90) sind(90) ; -sind(90) cosd(90)];
            oTPRProjData(curInds) = (targPlaneProj(curInds,:)-center_byRew(r,:))*(unitVec * rotMat)';
        end; clear r
    end; clear d

    % Z-score by direction if desired
    if contains(method,'dirZ')
        TPRProj_byDir = groupDataByLabel(TPRProjData(~skipLabels),[directionLabels(~skipLabels)]);
        ndirections = size(TPRProj_byDir,1); directions = unique(directionLabels);
        for d = 1:ndirections
            curInds = directionLabels==directions(d);
            TPRProjData(curInds) = (TPRProjData(curInds)-mean(TPRProj_byDir{d}))/std(TPRProj_byDir{d});
        end; clear d

        oTPRProj_byDir = groupDataByLabel(oTPRProjData(~skipLabels),[directionLabels(~skipLabels)]);
        for d = 1:ndirections
            curInds = directionLabels==directions(d);
            oTPRProjData(curInds) = (oTPRProjData(curInds)-mean(oTPRProj_byDir{d}))/std(oTPRProj_byDir{d});
        end; clear d
    end
    TPRProj_subjSplit{f} = TPRProjData(~skipLabels);
    oTPRProj_subjSplit{f} = oTPRProjData(~skipLabels);
    PSLs_subjSplit{f} = rotPSLocs(~skipLabels,:);
    directionLabels_subjSplit{f} = directionLabels(~skipLabels);
    rewardLabels_subjSplit{f} = rewardLabels(~skipLabels);
    successLabels_subjSplit{f} = trialStatusLabels(~skipLabels)==1;
%     successLabels_subjSplit{f} = ismember(trialStatusLabels,[1 -22 -34])
    dayLabels_subjSplit{f} = dayLabels(~skipLabels);
    prewLabels_subjSplit{f} = prevRewardLabels(~skipLabels); 
    
    % Get correlations with RT; no need for z-scores here since it's within condition
    TPRProj_byDir_byRew = groupDataByLabel(TPRProjData(~skipLabels),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    RTs_byDir_byRew = groupDataByLabel(RTs(~skipLabels),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    RTs_byDir_byRew_subjSplit{f} = RTs_byDir_byRew;
    [corr_RT_TPR_byDR_subjSplit{f},p_RT_TPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),RTs_byDir_byRew,TPRProj_byDir_byRew);
    corrP_RT_TPR_byDR_subjSplit{f} = cellfun(@(x,y) corr(x,y,'type','pearson'),RTs_byDir_byRew,TPRProj_byDir_byRew);
    PSs_byDir_byRew = groupDataByLabel(PSs(~skipLabels),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    [corr_PS_TPR_byDR_subjSplit{f},p_PS_TPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),PSs_byDir_byRew,TPRProj_byDir_byRew);
    HTs_byDir_byRew = groupDataByLabel(HTs(~skipLabels),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    [corr_HT_TPR_byDR_subjSplit{f},p_HT_TPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x(~isnan(x+y)),y(~isnan(x+y)),'type','spearman'),HTs_byDir_byRew,TPRProj_byDir_byRew);
    BEPs_byDir_byRew = groupDataByLabel(rotBREPs(~skipLabels,1),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    [corr_BEP_TPR_byDR_subjSplit{f},p_BEP_TPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),BEPs_byDir_byRew,TPRProj_byDir_byRew);
    absoffPSL_byDir_byRew = groupDataByLabel(abs(rotPSLocs(~skipLabels,2)),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    [corr_absoffPSL_TPR_byDR_subjSplit{f},p_absoffPSL_TPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),absoffPSL_byDir_byRew,TPRProj_byDir_byRew);
    
    % Same but for off-TPA
    oTPRProj_byDir_byRew = groupDataByLabel(oTPRProjData(~skipLabels),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    [corr_RT_oTPR_byDR_subjSplit{f},p_RT_oTPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),RTs_byDir_byRew,oTPRProj_byDir_byRew);
    corrP_RT_oTPR_byDR_subjSplit{f} = cellfun(@(x,y) corr(x,y,'type','pearson'),RTs_byDir_byRew,oTPRProj_byDir_byRew);
    [corr_PS_oTPR_byDR_subjSplit{f},p_PS_oTPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),PSs_byDir_byRew,oTPRProj_byDir_byRew);
    [corr_HT_oTPR_byDR_subjSplit{f},p_HT_oTPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x(~isnan(x+y)),y(~isnan(x+y)),'type','spearman'),HTs_byDir_byRew,oTPRProj_byDir_byRew);
    [corr_BEP_oTPR_byDR_subjSplit{f},p_BEP_oTPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),BEPs_byDir_byRew,oTPRProj_byDir_byRew);
    [corr_absoffPSL_oTPR_byDR_subjSplit{f},p_absoffPSL_oTPR_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),absoffPSL_byDir_byRew,oTPRProj_byDir_byRew);
    

    % Get correlation of PS (linearly conditioned on RT)
    RTs_byDir_byRew = groupDataByLabel(RTs(~skipLabels),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    betas_byDir_byRew = cellfun(@(x,y) mvregress([x ones(length(x),1)],y), RTs_byDir_byRew, PSs_byDir_byRew, 'uniformoutput', false);
    predPSs_byDir_byRew = cellfun(@(x,y) [x ones(length(x),1)]*y, RTs_byDir_byRew, betas_byDir_byRew,'uniformoutput', false);
    PSgRTs_byDir_byRew = cellfun(@(x,y) x-y, PSs_byDir_byRew, predPSs_byDir_byRew, 'uniformoutput', false);
    corr_PSgRT_TPR_byDR_subjSplit{f} = cellfun(@(x,y) corr(x,y,'type','spearman'),PSgRTs_byDir_byRew,TPRProj_byDir_byRew);
    
    % Do reverse as well
    betas_byDir_byRew = cellfun(@(x,y) mvregress([x ones(length(x),1)],y), PSs_byDir_byRew, RTs_byDir_byRew, 'uniformoutput', false);
    predRTs_byDir_byRew = cellfun(@(x,y) [x ones(length(x),1)]*y, PSs_byDir_byRew, betas_byDir_byRew,'uniformoutput', false);
    RTgPSs_byDir_byRew = cellfun(@(x,y) x-y, RTs_byDir_byRew, predRTs_byDir_byRew, 'uniformoutput', false);
    corr_RTgPS_TPR_byDR_subjSplit{f} = cellfun(@(x,y) corr(x,y,'type','spearman'),RTgPSs_byDir_byRew,TPRProj_byDir_byRew);
    
    disp('Corr between PS and RT')
    temp = cellfun(@(x,y) corr(x,y,'type','spearman'),RTs_byDir_byRew,PSs_byDir_byRew)
    mean(temp(:))



    % NEW THING: look at relationship between TPA/oTPA and endpt / PSL
    onPSLs_byDir_byRew = groupDataByLabel(rotPSLocs(~skipLabels,1),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    offPSLs_byDir_byRew = groupDataByLabel(rotPSLocs(~skipLabels,2),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    onEndpts_byDir_byRew = groupDataByLabel(rotEndpts(~skipLabels,1),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    offEndpts_byDir_byRew = groupDataByLabel(rotEndpts(~skipLabels,2),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);

    [corr_TPR_onPSL_byDR_subjSplit{f},p_TPR_onPSL_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),TPRProj_byDir_byRew,onPSLs_byDir_byRew);
    [corr_oTPR_offPSL_byDR_subjSplit{f},p_oTPR_offPSL_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman'),oTPRProj_byDir_byRew,offPSLs_byDir_byRew);
    [corr_TPR_onEndpt_byDR_subjSplit{f},p_TPR_onEndpt_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman','rows','complete'),TPRProj_byDir_byRew,onEndpts_byDir_byRew);
    [corr_oTPR_offEndpt_byDR_subjSplit{f},p_oTPR_offEndpt_byDR_subjSplit{f}] = cellfun(@(x,y) corr(x,y,'type','spearman','rows','complete'),oTPRProj_byDir_byRew,offEndpts_byDir_byRew);
    

end; clear f
disp('Done processing data')


%% NEW PLOTTING METHOD: combine over animals to make one histogram
rewColors = getDistinctColors('SELECT_ORDER',8);
nbins = 19;
% axisMax = ceil(max(abs(curFullData*10)))/10; % get outer tenth
axisMax = 0.75;
maxy = 30;
fs = 10;


figure; set(gcf,'position', [-586 29 293 954])

subplot(5,1,1); hold on
curFullData = cellArrayToVector(corr_RT_TPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_RT_TPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{RT,TPA}')
disp('RT')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
title('Corr. with TPA projections')
ylabel({'Reaction Time','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])

subplot(5,1,2); hold on
curFullData = cellArrayToVector(corr_PS_TPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_PS_TPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{PS,TPA}')
disp('PS')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Peak Speed','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])


subplot(5,1,3); hold on
curFullData = cellArrayToVector(corr_HT_TPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_HT_TPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{HT,TPA}')
disp('HT')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Homing Time','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])

subplot(5,1,4); hold on
curFullData = cellArrayToVector(corr_BEP_TPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_BEP_TPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{BEP,TPA}')
disp('BEP')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Ball. Endpt. Pred.','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])

subplot(5,1,5); hold on
curFullData = cellArrayToVector(corr_absoffPSL_TPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
pFullData = cellArrayToVector(p_absoffPSL_TPR_byDR_subjSplit);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
xlabel('Within-Condition Spearman Corr.')
disp('AE')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'|Off-Axis Error|','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])


% Save 
figname = ['FigS4_TPAcorrWithBehavior_histograms_' method]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);



%% SAME IDEA but for off-TPA


figure; set(gcf,'position', [-586 29 293 954])

subplot(5,1,1); hold on
curFullData = cellArrayToVector(corr_RT_oTPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_RT_oTPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{RT,TPA}')
disp('RT')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
title('Corr. with off-TPA projections')
ylabel({'Reaction Time','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])

subplot(5,1,2); hold on
curFullData = cellArrayToVector(corr_PS_oTPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_PS_oTPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{PS,TPA}')
disp('PS')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Peak Speed','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])


subplot(5,1,3); hold on
curFullData = cellArrayToVector(corr_HT_oTPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_HT_oTPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{HT,TPA}')
disp('HT')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Homing Time','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])

subplot(5,1,4); hold on
curFullData = cellArrayToVector(corr_BEP_oTPR_byDR_subjSplit);
pFullData = cellArrayToVector(p_BEP_oTPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
% xlabel('r_{BEP,TPA}')
disp('BEP')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Ball. Endpt. Pred.','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])

subplot(5,1,5); hold on
curFullData = cellArrayToVector(corr_absoffPSL_oTPR_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
pFullData = cellArrayToVector(p_absoffPSL_oTPR_byDR_subjSplit);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
xlabel('Within-Condition Spearman Corr.')
disp('AE')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
ylabel({'Off-Axis Error','Count'})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
axis([-axisMax axisMax 0 maxy])
xticks([-0.5 0 0.5])


% Save 
figname = ['FigS4_offTPAcorrWithBehavior_histograms_' method]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Make an example plot of the target plane versus peak speed location
f = 1;
d = 1;
r = 1;

curInds = rewardLabels_subjSplit{f}==r & directionLabels_subjSplit{f}==d;
curTPData = [TPRProj_subjSplit{f}(curInds) oTPRProj_subjSplit{f}(curInds)];
curPSLData = PSLs_subjSplit{f}(curInds,:);

figure; set(gcf,'position',[-1550 256 795 663])
ms = 10;
rewColors = getDistinctColors('SELECT_ORDER',8);
subplot(2,2,1); 
plot(curTPData(:,1), curTPData(:,2), '.','markersize',ms,'color',rewColors{r})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
xlabel('Target Preparation Axis (a.u.)')
ylabel('Off-Target Preparation Axis (a.u.)')
axis([-4 4 -4 4])

subplot(2,2,2); 
plot(curPSLData(:,1), curPSLData(:,2), '.','markersize',ms,'color',rewColors{r})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
xlabel('On-Axis Peak Speed Location (mm)')
ylabel('Off-Axis Peak Speed Location (mm)')
axis([20 45 -12.5 12.5])

subplot(2,2,3); 
plot(curPSLData(:,1), curTPData(:,1),'.','markersize',ms,'color',rewColors{r})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
xlabel('On-Axis Peak Speed Location (mm)')
ylabel('Target Preparation Axis (a.u.)')
axis([20 45 -4 4])
corr(curPSLData(:,1), curTPData(:,1),'type','spearman')

subplot(2,2,4); 
plot(curPSLData(:,2), curTPData(:,2),'.','markersize',ms,'color',rewColors{r})
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');
xlabel('Off-Axis Peak Speed Location (mm)')
ylabel('Off-Target Preparation Axis (a.u.)')
axis([-12.5 12.5 -4 4])
corr(curPSLData(:,2), curTPData(:,2),'type','spearman')

% Save 
figname = ['FigS4_TPandPSLScatters_' method]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);




%% Show correlation histograms

figure; set(gcf,'position',[-1213 543 624 403])
subplot(2,2,1);  hold on
curFullData = cellArrayToVector(corr_TPR_onPSL_byDR_subjSplit);
pFullData = cellArrayToVector(p_TPR_onPSL_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
ylabel('Count')
title('TPA v. On-Axis Peak Speed Location')
xlabel('Within-condition Spearman Corr.')
disp('On, PSL')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');

subplot(2,2,2);  hold on
curFullData = cellArrayToVector(corr_oTPR_offPSL_byDR_subjSplit);
pFullData = cellArrayToVector(p_oTPR_offPSL_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
title('Off-TPA v. Off-Axis Peak Speed Location')
xlabel('Within-condition Spearman Corr.')
disp('Off, PSL')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');

subplot(2,2,3);  hold on
curFullData = cellArrayToVector(corr_TPR_onEndpt_byDR_subjSplit);
pFullData = cellArrayToVector(p_TPR_onEndpt_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
ylabel('Count')
title('TPA v. On-Axis Reach Endpoint')
xlabel('Within-condition Spearman Corr.')
disp('On, endpt')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');

subplot(2,2,4);  hold on
curFullData = cellArrayToVector(corr_oTPR_offEndpt_byDR_subjSplit);
pFullData = cellArrayToVector(p_oTPR_offEndpt_byDR_subjSplit);
binEdges = linspace(-axisMax,axisMax,nbins);
histogram(curFullData,binEdges,'facecolor',[0.85 0.85 0.85]);
[vals,xvals] = histcounts(curFullData,binEdges);
curSigData = curFullData(pFullData < 0.05);
histogram(curSigData,binEdges,'facecolor',[0.35 0.35 0.35]);
plot(median(curFullData)*ones(100,1),linspace(0,maxy),'k--','linewidth',1.5)
plot(zeros(100,1),linspace(0,maxy),'k-','linewidth',1.5)
title('Off-TPA v. Off-Axis Reach Endpoint')
xlabel('Within-condition Spearman Corr.')
disp('Off, endpt')
disp(['median = ' num2str(median(curFullData))])
disp(['p = ' num2str(signrank(curFullData))])
set(gca,'fontsize',fs,'fontname','arial','tickdir','out');



% Save 
figname = ['FigS4_TPandAccHistograms_' method]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

