%%% This script makes figures related to EMG and reward, namely:
%%% - average EMG traces for the whole trial as a function of reward
%%% - the delay EMG directional tuning curve as a function of reward
%%% - comparisons between delay EMG and the reward axis
%%%
%%% Adam Smoulder, 1/18/23

neuralAndBehFilename = 'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
EMGFilename = 'RockyPreprocEMGData_normalized.mat';
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;
rng(3195)

%% Load all data
load(neuralAndBehFilename)
load(EMGFilename);

%% Get average EMG traces for each trial
timeTO = -50:250; ntimeTO = length(timeTO); % Aligned to target onset
% timeMO = -200:200; ntimeMO = length(timeMO); % aligned to movement onset
timeMO = -300:500; ntimeMO = length(timeMO); % aligned to GC
timeTA = -50:250; ntimeTA = length(timeTA); % aligned to target acquisition
kernelSize = 51; % size of smoothing kernel
kernel = ones(kernelSize,1)/kernelSize;
HKW = floor(kernelSize/2); % half kernel width = used for padding edges with extra data

validMuscleLabels = [EMGData.validMuscleLabel]';
[ntrials,nmuscles] = size(validMuscleLabels);
EMG_TO_byTrial_byTime_byMuscle = nan(ntrials,ntimeTO,nmuscles);
EMG_MO_byTrial_byTime_byMuscle = nan(ntrials,ntimeMO,nmuscles);
EMG_TA_byTrial_byTime_byMuscle = nan(ntrials,ntimeTA,nmuscles);
for i = 1:ntrials
    if trialStatusLabels(i)==1 && delayLengths(i) >= 317 % Rocky's minimum delay for choking stuff
        % Get event times
        stateTable = EMGData(i).stateTable;
        TO = stateTable(2,stateTable(1,:)==3); % state 3 is target onset
%         MO = stateTable(2,stateTable(1,:)==4)+RTs(i); % state 4 is go cue
        MO = stateTable(2,stateTable(1,:)==4); % state 4 is go cue
        TA = stateTable(2,stateTable(1,:)==6); % state 6 is target acq or failure
        
        % Get the EMG window around each of these, smooth, then store
        TOInds = ismember(EMGData(i).time,TO+(min(timeTO)-HKW:max(timeTO)+HKW));
        EMG_TO_byTrial_byTime_byMuscle(i,:,:) = convByDim((EMGData(i).EMG(TOInds,:)),kernel,1,'valid');
        MOInds = ismember(EMGData(i).time,MO+(min(timeMO)-HKW:max(timeMO)+HKW));
        EMG_MO_byTrial_byTime_byMuscle(i,:,:) = convByDim((EMGData(i).EMG(MOInds,:)),kernel,1,'valid');
        TAInds = ismember(EMGData(i).time,TA+(min(timeTA)-HKW:max(timeTA)+HKW));
        EMG_TA_byTrial_byTime_byMuscle(i,:,:) = convByDim((EMGData(i).EMG(TAInds,:)),kernel,1,'valid');
    end
end; clear i

% Go through each muscle and replace bad day indices with nans
for m = 1:nmuscles
    EMG_TO_byTrial_byTime_byMuscle(~validMuscleLabels(:,m),:,m) = nan;
    EMG_MO_byTrial_byTime_byMuscle(~validMuscleLabels(:,m),:,m) = nan;
    EMG_TA_byTrial_byTime_byMuscle(~validMuscleLabels(:,m),:,m) = nan;
end; clear m

%% Get average EMG before target onset and around go cue
TOTime = -200:0;
GCTime = -150:50;
EMG_pTO = nan(ntrials,nmuscles);
EMG_pGC = nan(ntrials,nmuscles);
for i = 1:ntrials
    if delayLengths(i) >= 317 % Rocky's minimum delay for choking stuff
        % Get event times
        stateTable = EMGData(i).stateTable;
        TO = stateTable(2,stateTable(1,:)==3); % state 3 is target onset
        EMG_pTO(i,:) = mean(EMGData(i).EMG(ismember(EMGData(i).time, TO+TOTime),:));
        GC = stateTable(2,stateTable(1,:)==4); % state 4 is go cue
        EMG_pGC(i,:) = mean(EMGData(i).EMG(ismember(EMGData(i).time, GC+GCTime),:));
    end
end; clear i


%% Make a plot showing the average EMG trace by reward for an example direction
d = 8;
clear EMG_TO_dirMean_byTime_byMuscle_byRew EMG_TO_dirMedian_byTime_byMuscle_byRew EMG_TO_dirSem_byTime_byMuscle_byRew
clear EMG_MO_dirMean_byTime_byMuscle_byRew EMG_MO_dirMedian_byTime_byMuscle_byRew EMG_MO_dirSem_byTime_byMuscle_byRew
clear EMG_TA_dirMean_byTime_byMuscle_byRew EMG_TA_dirMedian_byTime_byMuscle_byRew EMG_TA_dirSem_byTime_byMuscle_byRew
for r = 1:nrewards
    curInds = directionLabels==d & rewardLabels==r;
    
    EMG_TO_dirMean_byTime_byMuscle_byRew(:,:,r) = squeeze(nanmean(EMG_TO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TO_dirMedian_byTime_byMuscle_byRew(:,:,r) = squeeze(nanmedian(EMG_TO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TO_dirSem_byTime_byMuscle_byRew(:,:,r) = squeeze(nansem(EMG_TO_byTrial_byTime_byMuscle(curInds,:,:),1));
    
    EMG_MO_dirMean_byTime_byMuscle_byRew(:,:,r) = squeeze(nanmean(EMG_MO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_MO_dirMedian_byTime_byMuscle_byRew(:,:,r) = squeeze(nanmedian(EMG_MO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_MO_dirSem_byTime_byMuscle_byRew(:,:,r) = squeeze(nansem(EMG_MO_byTrial_byTime_byMuscle(curInds,:,:),1));
    
    EMG_TA_dirMean_byTime_byMuscle_byRew(:,:,r) = squeeze(nanmean(EMG_TA_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TA_dirMedian_byTime_byMuscle_byRew(:,:,r) = squeeze(nanmedian(EMG_TA_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TA_dirSem_byTime_byMuscle_byRew(:,:,r) = squeeze(nansem(EMG_TA_byTrial_byTime_byMuscle(curInds,:,:),1));
end; clear r

% %%
% gap = nan(49,nmuscles); % gap between traces
% fullMeanTrace = [EMG_TO_dirMean_byTime_byMuscle_byRew(:,:,1) ; gap ; 
%                  EMG_MO_dirMean_byTime_byMuscle_byRew(:,:,1) ; gap ; 
%                  EMG_TA_dirMean_byTime_byMuscle_byRew(:,:,1)];
% ntimefull = length(fullMeanTrace);
% figure; plot(fullMeanTrace);
% xticks(1:50:ntimefull)
% xticklabels([min(timeTO):50:max(timeTO) min(timeMO):50:max(timeMO) min(timeTA):50:max(timeTA)])

% Make a figure
gap = nan(49,1); % gap between traces
concatTime = [min(timeTO):50:max(timeTO) min(timeMO):50:max(timeMO) min(timeTA):50:max(timeTA)];
zeroInds = find(concatTime==0);
figure
rewColors = getDistinctColors('SELECT_ORDER',8);
lw = 1.5;
for m = 1:nmuscles
    subplot(1,nmuscles,m); hold on
    for r = 1:nrewards
        fullMeanTrace = [EMG_TO_dirMean_byTime_byMuscle_byRew(:,m,r) ; gap ;
            EMG_MO_dirMean_byTime_byMuscle_byRew(:,m,r) ; gap ;
            EMG_TA_dirMean_byTime_byMuscle_byRew(:,m,r)] - nanmean(EMG_pTO(:,m)); % a.u., so set each muscle 0 to pTO
        ntimefull = length(fullMeanTrace);
        plot(fullMeanTrace,'-','color',rewColors{r},'linewidth',lw)
        indVals = 1:50:ntimefull;
        xticks(indVals(zeroInds))
%         xticklabels({'TO','MO','TA'})
        xticklabels({'TO','GC','TA'}) % use this if we align to GC instead
%         xticks()
%         xticklabels([min(timeTO):50:max(timeTO) min(timeMO):50:max(timeMO) min(timeTA):50:max(timeTA)])
        
    end; clear r
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    title(muscleNames{m})
    if m == 1
        ylabel(num2str(d))
    end
end; clear m
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 1924 338])

%% Same idea, but don't care about reward, just show directions
clear EMG_TO_dirMean_byTime_byMuscle_byDir EMG_TO_dirMedian_byTime_byMuscle_byDir EMG_TO_dirSem_byTime_byMuscle_byDir
clear EMG_MO_dirMean_byTime_byMuscle_byDir EMG_MO_dirMedian_byTime_byMuscle_byDir EMG_MO_dirSem_byTime_byMuscle_byDir
clear EMG_TA_dirMean_byTime_byMuscle_byDir EMG_TA_dirMedian_byTime_byMuscle_byDir EMG_TA_dirSem_byTime_byMuscle_byDir
for d = 1:ndirections
    curInds = directionLabels==d;
    
    EMG_TO_dirMean_byTime_byMuscle_byDir(:,:,d) = squeeze(nanmean(EMG_TO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TO_dirMedian_byTime_byMuscle_byDir(:,:,d) = squeeze(nanmedian(EMG_TO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TO_dirSem_byTime_byMuscle_byDir(:,:,d) = squeeze(nansem(EMG_TO_byTrial_byTime_byMuscle(curInds,:,:),1));
    
    EMG_MO_dirMean_byTime_byMuscle_byDir(:,:,d) = squeeze(nanmean(EMG_MO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_MO_dirMedian_byTime_byMuscle_byDir(:,:,d) = squeeze(nanmedian(EMG_MO_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_MO_dirSem_byTime_byMuscle_byDir(:,:,d) = squeeze(nansem(EMG_MO_byTrial_byTime_byMuscle(curInds,:,:),1));
    
    EMG_TA_dirMean_byTime_byMuscle_byDir(:,:,d) = squeeze(nanmean(EMG_TA_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TA_dirMedian_byTime_byMuscle_byDir(:,:,d) = squeeze(nanmedian(EMG_TA_byTrial_byTime_byMuscle(curInds,:,:)));
    EMG_TA_dirSem_byTime_byMuscle_byDir(:,:,d) = squeeze(nansem(EMG_TA_byTrial_byTime_byMuscle(curInds,:,:),1));
end; clear d


gap = nan(49,1); % gap between traces
concatTime = [min(timeTO):50:max(timeTO) min(timeMO):50:max(timeMO) min(timeTA):50:max(timeTA)];
zeroInds = find(concatTime==0);


figure
set(gcf,'renderer','painters')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[0 500] 324 954])
dirColors = getDistinctColors('SELECT_ORDER',16);
lw = 1.5;
for m = 1:nmuscles
    subplot(nmuscles,1,m); hold on
    for d = 1:ndirections
        y = [EMG_TO_dirMean_byTime_byMuscle_byDir(:,m,d) ; gap ;
            EMG_MO_dirMean_byTime_byMuscle_byDir(:,m,d) ; gap ;
            EMG_TA_dirMean_byTime_byMuscle_byDir(:,m,d)] - nanmean(EMG_pTO(:,m)); % a.u., so set each muscle 0 to pTO
        err = [EMG_TO_dirSem_byTime_byMuscle_byDir(:,m,d) ; gap ;
            EMG_MO_dirSem_byTime_byMuscle_byDir(:,m,d) ; gap ;
            EMG_TA_dirSem_byTime_byMuscle_byDir(:,m,d)]; % a.u., so set each muscle 0 to pTO
        ntimefull = length(fullMeanTrace);
%         plot(fullMeanTrace,'-','color',dirColors{d},'linewidth',lw)
        x = 1:ntimefull;
        shadedErrorBar(x,y,err,...
           'lineProps',{'r-','color',dirColors{d},'linewidth',lw})
        indVals = 1:50:ntimefull;
        xticks(indVals(zeroInds))
        %         xticklabels({'TO','MO','TA'})
        %         xticklabels({'TO','GC','TA'}) % use this if we align to GC instead
        xticklabels([])
        
    end; clear r
    set(gca,'fontsize',10,'fontname','arial','tickdir','out','box','off');
    %     title(muscleNames{m})
    %     axis([-inf inf -inf inf])
%     axis([-inf inf -3 8])
    axis([-inf inf -inf inf])
    %     if m == 1
    ylabel('EMG (a.u.)')
    %     end
    if m == nmuscles
        xlabel(['Time in trial'])
    end
end; clear m

% Save it!
figname = ['FigS1_EMGTracesByDirection_errbar' ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

%% Plot delay period EMG tuning curves by reward

% Evaluate if there is any reward tuning for each muscle using 2-way anova
pvals_byCond_byMuscle = nan(3,nmuscles); % dir, rew, interaction
EMGMean_byDir_byRew_byMuscle = nan(ndirections,nrewards,nmuscles);
EMGSem_byDir_byRew_byMuscle = nan(ndirections,nrewards,nmuscles);
for m = 1:nmuscles
    X = EMG_pGC(:,m)-EMG_pTO(:,m);
    Y1 = directionLabels;
    Y2 = rewardLabels;
    badInds = isnan(X);
    Y1(badInds) = []; Y2(badInds) = []; X(badInds) = [];
    [pvals_byCond_byMuscle(:,m),~,stats] = anovan(X,[Y1 Y2],'display','off','model','interaction');
    
    % Get the mean and SEM for plots
    curEMG_byDir_byRew = groupDataByLabel(X,[Y1 Y2]);
    EMGMean_byDir_byRew_byMuscle(:,:,m) = cellfun(@(x) mean(x), curEMG_byDir_byRew);
    EMGSem_byDir_byRew_byMuscle(:,:,m) = cellfun(@(x) nansem(x,1), curEMG_byDir_byRew);
end; clear m

% Make a figure showing these
figure
rewColors = getDistinctColors('SELECT_ORDER',8);
for m = 1:nmuscles
    subplot(nmuscles,1,m); hold on
    for r = 1:nrewards
        errorbar((1:ndirections)-0.25+0.1*r, EMGMean_byDir_byRew_byMuscle(:,r,m), EMGSem_byDir_byRew_byMuscle(:,r,m), ...
            '-','linewidth',lw,'color',rewColors{r},'capsize',1)
    end; clear r
    set(gca,'fontsize',10,'fontname','arial','tickdir','out','box','off');
    %     title(muscleNames{m})
    axis([0.5 0.5+ndirections -inf inf])
    %     if m == 1
    ylabel('EMG (a.u.)')
    %     end
    if m == 5
        xlabel(['Direction'])
    end
    xticks(1:ndirections)
    if m == 5
        xticklabels({['0' char(176)],['45' char(176)],['90' char(176)],['135' char(176)],['180' char(176)],['225' char(176)],['270' char(176)],['315' char(176)]})
    else
        xticklabels([])
    end
end; clear m
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[0 500] 324 954])

% Save it!
figname = ['FigS1_delayEMGByReward' ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);