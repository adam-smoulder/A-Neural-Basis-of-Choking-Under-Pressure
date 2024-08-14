%%% This figure looks at the behavior of the animals as a function of
%%% reward.
%%% 
%%% RT: U for E, almost U for P, U with high J for R. Var is U for all, for
%%% Prez J is very high
%%% PS: mean inv-U, but to varying degree (prez is low J, rocky is mostly low
%%% J, Earl is symmetric). Variance is Small-high
%%% HT: Median is mostly mono-up, var is inconsistent
%%% BREP: mean is strong mono down for E, but U for R/P. Var is U for all
%%%
%%% Adam Smoulder, 10/26/23 (edit 4/17/24)


filename_bySubject = {...
    'Fig1BehaviorData_Earl_20240522_125841';
    'Fig1BehaviorData_Prez_20240522_125958';
    'Fig1BehaviorData_Rocky_20240522_130140';
    ...'Fig1BehaviorData_Nelson_20221208_131310.mat';
    ...'Fig1BehaviorData_Ford_20221208_131312.mat';
    };
nsubjects = length(filename_bySubject);
figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;


%% For each animal, we'll get behavioral metrics
nboots = 10000;
epochGroups = {[1 -11 -12 -22 -23 -34],[1 -22 -23 -34],[1 -34]}; % for consideration in delay, reach, and hold
epochSuccessGroups = {[1 -22 -23 -34],[1 -34],[1]};
targetDistance_bySubj = [85 85 85 65 80]; % diff for each subject
nepochs = length(epochGroups);

statusRate_subjSplit_epochSplit_statusSplit_byRew = cell(nsubjects,nepochs,3); % 3 is the max # of unique statuses in one epoch
pvals_subjSplit_epochSplit_statusSplit_byRewXRew = cell(nsubjects,nepochs,3); 
succRate_subjSplit_byRew_byDiff = cell(nsubjects,1);
succRate_subjSplit_byRew_byDiff_sem = cell(nsubjects,1);
succRate_subjSplit_byRew_byPSF = cell(nsubjects,1);
succRate_subjSplit_byRew_byPSF_sem = cell(nsubjects,1);
succRate_subjSplit_byPrew = cell(nsubjects,1);
succRate_subjSplit_byPrew_sem = cell(nsubjects,1);
RT_subjSplit_byRew = cell(nsubjects,1);
PS_subjSplit_byRew = cell(nsubjects,1);
HT_subjSplit_byRew = cell(nsubjects,1);
BREP_subjSplit_byRew = cell(nsubjects,1);
% endPt_subjSplit_byRew = cell(nsubjects,1);
% endptErr_subjSplit_byRew = cell(nsubjects,1);
dirLabels_subjSplit_metricSplit_byRew = cell(nsubjects,6);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = subjectName;
    badInds = false(length(directionLabels),1); tag = [];
%     badInds = delayLengths < (400-8); tag = '_noShortDelay'% skip short delay inds
    ntrials = sum(~badInds);

    % Status rate by epoch
    for e = 1:nepochs
        % Get the number of trials for each status in this epoch
        curInds = ismember(trialStatusLabels,epochGroups{e}) & ~badInds;
        epochStatusLabels = trialStatusLabels;
        epochStatusLabels(ismember(epochStatusLabels,epochSuccessGroups{e})) = 1; % it's a success for this epoch
        n_byRew_byStatus = nByLabel([rewardLabels(curInds) epochStatusLabels(curInds)]);
        
        % Get rate for each status and p-values
        curStatuses = unique(epochStatusLabels(curInds));
        for s = 1:length(curStatuses)
            statusRate_subjSplit_epochSplit_statusSplit_byRew{f,e,s} = 100*n_byRew_byStatus(:,s)./sum(n_byRew_byStatus,2);
            [pvals_subjSplit_epochSplit_statusSplit_byRewXRew{f,e,s}] = binomialProportionTest2(n_byRew_byStatus(:,s)./sum(n_byRew_byStatus,2),sum(n_byRew_byStatus,2));
        end; clear s
    end; clear e
    
    % Behavioral metrics by reward (get summary stats later)
    validRTInds = ~badInds & ~isnan(RTs) & RTs >= 100 & RTs <= 550;
    RT_subjSplit_byRew{f} = groupDataByLabel(RTs(validRTInds),rewardLabels(validRTInds));
    validPSInds = ~badInds & ~isnan(PSs) & PSs >= 0.1;
    PS_subjSplit_byRew{f} = groupDataByLabel(PSs(validPSInds),rewardLabels(validPSInds));
    validHTInds = ~badInds & ~isnan(HTs);
    HT_subjSplit_byRew{f} = groupDataByLabel(HTs(validHTInds),rewardLabels(validHTInds));
    BREP_subjSplit_byRew{f} = groupDataByLabel(rotBREPs(validPSInds,1),rewardLabels(validPSInds)); % We'll just use on-axis to start
    validAngleInds = validPSInds & ~isnan(angleError_endpt_signed);
    % AA_subjSplit_byRew{f} = groupDataByLabel(abs(angleError_endpt_signed(validAngleInds)),rewardLabels(validAngleInds));
    % AA_subjSplit_byRew{f} = groupDataByLabel(abs(angleError_pred_signed(validAngleInds)),rewardLabels(validAngleInds)); % pred is based off of start and peak speed location; more "ballistic"-y
    offPSL_subjSplit_byRew{f} = groupDataByLabel(rotPSLocs(validAngleInds,2),[rewardLabels(validAngleInds)]);

    % Also get direction labels for each of these
    dirLabels_subjSplit_metricSplit_byRew{f,1} = groupDataByLabel(directionLabels(validRTInds),rewardLabels(validRTInds));
    dirLabels_subjSplit_metricSplit_byRew{f,2} = groupDataByLabel(directionLabels(validPSInds),rewardLabels(validPSInds));
    dirLabels_subjSplit_metricSplit_byRew{f,3} = groupDataByLabel(directionLabels(validHTInds),rewardLabels(validHTInds));
    dirLabels_subjSplit_metricSplit_byRew{f,4} = groupDataByLabel(directionLabels(validPSInds),rewardLabels(validPSInds));
    dirLabels_subjSplit_metricSplit_byRew{f,5} = groupDataByLabel(directionLabels(validAngleInds),rewardLabels(validAngleInds));
end; clear f


%% Make plots for the status rates by epoch
statusColors_epochSplit = {{[0 0.85 0],[1 0.3 0],[1 0.65 0]},{[0 0.85 0],[0.3 0.3 0.3],[0.8 0.8 0.8]},{[0 0.85 0],[0.9 0.4 0.9]}};
ymin_byEpoch = [65 45 80];
lw = 1.5;
rewNames = {'S','M','L','J'};
figure
for f = 1:nsubjects
    for e = 1:nepochs
        subplot(nepochs,nsubjects,f+(e-1)*nsubjects); hold on
        curData = fliplr([statusRate_subjSplit_epochSplit_statusSplit_byRew{f,e,:}]); % rew x status (succ is first col)
        curColors = statusColors_epochSplit{e};
        if e ~= 3 % rearrange for visual appeal
            curData = curData(:,[1 3 2]);
            curColors = curColors([1 3 2]);
        end
        p = bar(curData,'stacked','linewidth',lw);
        axis([0 nrewards+1 ymin_byEpoch(e) 100])
        for s = 1:size(curData,2)
            p(s).FaceColor = curColors{s};
        end; clear s
        plot(curData(:,1),'k.-','linewidth',2,'markersize',20)
        
        % Label stuff
        set(gca,'fontsize',10,'fontname','arial','tickdir','out');
        xticks(1:nrewards); xticklabels(rewNames)
        if e == nepochs && f == ceil(nsubjects/2)
            xlabel('Reward')
        end
        if f == 1
            if e == 1
                ylabel('Delay Status Rate (%)')
            elseif e == 2
                ylabel('Reach Status Rate (%)')
            elseif e == 3
                ylabel('Target Hold Status Rate (%)')
            end
        else
            yticklabels([])
        end
        if e == 1
            title(['Monkey ' subjectNames{f}(1)])
        end
    end; clear e
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[0 550] 810 944])

% Save the figure
figname = ['FigS1_behaviorByReward_statusRatesByEpoch' tag];
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Show behavioral metrics
metricTags = {'RT','PS','HT','BREP','OAE'};
metricNames = {'Reaction Time (ms)','Peak Speed (m/s)','Homing Time (m/s)','Ball. Endpt. Pred (mm)','Off-Axis Deviation (mm)'}; %'Endpoint (mm)','Endpt error (mm)'};
metricColors = getDistinctColors('SELECT_ORDER',0);

figure
lw = 1.5;
ms = 20;

for m = 1:4 % show 3 per fig
    figure
    for f = 1:nsubjects
        % Get mean or median and var +/- SE
        eval(['y = ' metricTags{m} '_subjSplit_byRew{f};'])
        if any(cellfun(@(x) strcmp(metricTags{m},x),{'RT','HT','AA','endptErr'})) % not normally distributed; compare medians
            vals = cellfun(@(x) median(x),y);
            vals_err = cellfun(@(x) nansemed(x,1),y);
        else % use means
            vals = cellfun(@(x) mean(x),y);
            vals_err = cellfun(@(x) nansem(x,1),y);
        end
        
        subplot(1,5,f)
        errorbar(rewards,vals,vals_err,'k.-','color',metricColors{m},'linewidth',lw,'markersize',ms)
        axis([0.5 nrewards+0.5 -inf inf])
        
        % Label stuff
%         if m == 1
            title(['Monkey ' subjectNames{f}(1)])
%         end
        if f == 1
            ylabel(metricNames{m})
        end
        xticks(rewards); xticklabels(rewNames)
        if f == ceil(nsubjects/2)
            xlabel('Reward')
        end
        
        set(gca,'fontsize',10,'fontname','arial','tickdir','out','box','off');
    end; clear f
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1:2) 829 244])
    
    % Save the figure
    figname = ['FigS1_behaviorByReward_' metricTags{m} tag];
    saveFigAndSvg([figfolder 'recentFigs\'],figname);
    saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
end

%% For "angular error" we'll use the standard deviation of the peak speed
%  location off-axis error distributions.

nboots = 10000; % for error bars
AE_mean_subjSplit = cell(nsubjects,1);
AE_bootMean_subjSplit = cell(nsubjects,1);
AE_err_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    for r = 1:nrewards
        dirLabels = dirLabels_subjSplit_metricSplit_byRew{f,5}{r};
        offPSLs = offPSL_subjSplit_byRew{f}{r};
        [curOffPSLs_byDir,inds_byDir] = groupDataByLabel(offPSLs, dirLabels);
        ndirections = length(inds_byDir);
        AE_mean_subjSplit{f}(r) = mean(cellfun(@(x) std(x), curOffPSLs_byDir));
        AE_bootMeans_byDir = nan(nboots,ndirections);
        for b = 1:nboots
            bootInds = cellfun(@(x) x(randi(length(x),[length(x) 1])), inds_byDir, 'uniformoutput', false);
            bootOffPSLs_byDir = cellfun(@(x) offPSLs(x), bootInds, 'uniformoutput', false);
            AE_bootMeans_byDir(b,:) = cellfun(@(x) std(x), bootOffPSLs_byDir);
        end; clear b
        AE_bootMean_subjSplit{f}(r) = mean(AE_bootMeans_byDir(:));
        AE_err_subjSplit{f}(r) = mean(std(AE_bootMeans_byDir));
    end; clear r
    f
end; clear f


%% Make a figure
figure
m = 5;
for f = 1:nsubjects
    vals = AE_bootMean_subjSplit{f};
    err = AE_err_subjSplit{f};

    subplot(1,5,f)
    errorbar(rewards,vals,vals_err,'k.-','color',metricColors{m},'linewidth',lw,'markersize',ms)
    axis([0.5 nrewards+0.5 -inf inf])

    % Label stuff
    %         if m == 1
    title(['Monkey ' subjectNames{f}(1)])
    %         end
    if f == 1
        ylabel(metricNames{m})
    end
    xticks(rewards); xticklabels(rewNames)
    if f == ceil(nsubjects/2)
        xlabel('Reward')
    end

    set(gca,'fontsize',10,'fontname','arial','tickdir','out','box','off');
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 829 244])

% Save the figure
figname = ['FigS1_behaviorByReward_offAxisDeviation' tag];
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% If you want to do it with IQR instead, you can do this.

nboots = 1000; % for error bars
AE_mean_subjSplit = cell(nsubjects,1);
AE_bootMean_subjSplit = cell(nsubjects,1);
AE_err_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    for r = 1:nrewards
        dirLabels = dirLabels_subjSplit_metricSplit_byRew{f,5}{r};
        offPSLs = offPSL_subjSplit_byRew{f}{r};
        [curOffPSLs_byDir,inds_byDir] = groupDataByLabel(offPSLs, dirLabels);
        ndirections = length(inds_byDir);
        AE_mean_subjSplit{f}(r) = mean(cellfun(@(x) iqr(x), curOffPSLs_byDir));
        AE_bootMeans_byDir = nan(nboots,ndirections);
        for b = 1:nboots
            bootInds = cellfun(@(x) x(randi(length(x),[length(x) 1])), inds_byDir, 'uniformoutput', false);
            bootOffPSLs_byDir = cellfun(@(x) offPSLs(x), bootInds, 'uniformoutput', false);
            AE_bootMeans_byDir(b,:) = cellfun(@(x) iqr(x), bootOffPSLs_byDir);
        end; clear b
        AE_bootMean_subjSplit{f}(r) = mean(AE_bootMeans_byDir(:));
        AE_err_subjSplit{f}(r) = mean(std(AE_bootMeans_byDir));
    end; clear r
    f
end; clear f
