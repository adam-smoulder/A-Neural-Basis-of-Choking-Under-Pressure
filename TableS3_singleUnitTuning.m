%%% This script looks at single neuron tuning as a function of reward at
%%% the end of the delay period. We'll also look at directional tuning as
%%% well.
%%% 
%%% Adam Smoulder, 8/2/22 (edit 12/19/22)


addpath(genpath('D:\AdamMatlab\'))
filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
dateString = grabDateTimeString;
rng(3195)

%% For each animal, assess tuning and store statistics
% nmin_bySubj = [15 15 15]; % minimum number of trials per condition for analysis
nmin_bySubj = [10 10 10]; % minimum number of trials per condition for analysis

subjectNames = cell(nsubjects,1);
fracDataLost_bySubj = nan(nsubjects,1); % Amount of data lost if only keeping single neurons here
meanUnitFR_byDir_byRew_subjSplit = cell(nsubjects,1);
semUnitFR_byDir_byRew_subjSplit = cell(nsubjects,1);
varUnitFR_byDir_byRew_subjSplit = cell(nsubjects,1);
pval_Dir_Rew_subjSplit = cell(nsubjects,1);
rewTuningShapeLabel_subjSplit = cell(nsubjects,1);
pvalDir_byRew_subjSplit = cell(nsubjects,1);
pvalRew_byDir_subjSplit = cell(nsubjects,1);
pvalRew_byDir_rewComp_subjSplit = cell(nsubjects,1);
rewTuningShapeLabel_byDir_subjSplit = cell(nsubjects,1);
rewTuningChanges_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    [ntrials,~] = size(singleNeuronData_pGC);
    
    if strcmp(taskInfo.subjectName,'Earl')
        skipLabels = sum(isoutlier(binnedData_pGC,'thresholdfactor',4),2)>0; % VERY outlier-y, but present enough to influence...
    elseif strcmp(taskInfo.subjectName,'Rocky')
        skipLabels = delayLengths < 317;
    else % Prez's neural activity asymptotes in time, so no bad inds
        skipLabels = false(length(binnedData_pGC),1);
    end
    disp(['n outliers/to skip = ' num2str(sum(skipLabels))])
    singleNeuronData_pGC(skipLabels,:) = [];
    rewardLabels(skipLabels) = [];
    directionLabels(skipLabels) = [];
    
    % Identify units with at least nmin of each direction x reward
    % condition and throw away all else. 
    [validUnitInds,nminByUnit] = unitsThatHaveMinTrialCount(singleNeuronData_pGC,...
        nmin_bySubj(f),[directionLabels rewardLabels]);
    nunits = sum(validUnitInds);
    validNeuronData = singleNeuronData_pGC(:,validUnitInds);
    
    % Calculate how much data is thrown out with this method
    fracDataLost_bySubj(f) = 1-sum(~isnan(validNeuronData(:)))/sum(~isnan(singleNeuronData_pGC(:)));
    
    % For the remaining valid units, get tuning properties
    for u = 1:nunits
        curData = validNeuronData(:,u);
        validInds = ~isnan(curData);
        meanUnitFR_byDir_byRew_subjSplit{f}(u,:,:) = cellfun(@(x) mean(x), groupDataByLabel(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)]));
        semUnitFR_byDir_byRew_subjSplit{f}(u,:,:) = cellfun(@(x) nansem(x,1), groupDataByLabel(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)]));
        varUnitFR_byDir_byRew_subjSplit{f}(u,:,:) = cellfun(@(x) var(x), groupDataByLabel(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)]));
        
%         % Use 2-way ANOVA to assess if there is any dir or rew tuning
% %         pval_Dir_Rew_subjSplit{f}(u,:) = anovan(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)],'display','off');
%         [pval_Dir_Rew_subjSplit{f}(u,:),~,stats] = anovan(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)],'display','off','model','interaction');
%         statCompResults = multcompare(stats,'dim',2,'display','off');
%         pwRewPvals = statCompResults(:,end);
%         rewMeans = squeeze(mean(meanUnitFR_byDir_byRew_subjSplit{f}(u,:,:),2));
%         [shapeLabel,changes,shapeNames] = evalRewardTuningShape(pwRewPvals,rewMeans);
        
        % Trying with combining M/L for simplicity
        tempRewardLabels = rewardLabels;
        tempRewardLabels(rewardLabels==3) = 2;
        tempRewardLabels(rewardLabels==4) = 3;
        [pval_Dir_Rew_subjSplit{f}(u,:),~,stats] = anovan(curData(validInds),[directionLabels(validInds) tempRewardLabels(validInds)],'display','off','model','interaction');
        statCompResults = multcompare(stats,'dim',2,'display','off');
        pwRewPvals = statCompResults(:,end);
        rewMeans = mean(cellfun(@(x) mean(x), groupDataByLabel(curData(validInds),[directionLabels(validInds) tempRewardLabels(validInds)])))';
        [shapeLabel,changes,shapeNames] = evalRewardTuningShape3Rew(pwRewPvals,rewMeans);
        rewTuningShapeLabel_subjSplit{f}(u) = shapeLabel;
        
        % For each reward, assess direction tuning with KW test
        for r = 1:nrewards
            curInds = rewardLabels==rewards(r) & validInds;
            [pvalDir_byRew_subjSplit{f}(u,r),~,stats] = kruskalwallis(curData(curInds),directionLabels(curInds),'off');
        end; clear r
        
        % For each direction, assess reward tuning with KW test
        for d = 1:ndirections
            curInds = directionLabels==directions(d) & validInds;
            [pvalRew_byDir_subjSplit{f}(u,d),~,stats] = kruskalwallis(curData(curInds),rewardLabels(curInds),'off');
            T = multcompare(stats,'display','off');
            pvalRew_byDir_rewComp_subjSplit{f}(u,d,:) = T(:,end);
        end; clear d
        
        % Assess the tuning shape for each of these
        for d = 1:ndirections
            if pvalRew_byDir_subjSplit{f}(u,d) < 0.05
                p = squeeze(pvalRew_byDir_rewComp_subjSplit{f}(u,d,:));
                rewMeans = squeeze(meanUnitFR_byDir_byRew_subjSplit{f}(u,d,:));
                [shapeLabel,changes,shapeNames] = evalRewardTuningShape(p,rewMeans);
                rewTuningChanges_subjSplit{f}(u,d,:) = changes;
            else % no sig change; force 'None')
                shapeLabel = 6; % none
            end
            rewTuningShapeLabel_byDir_subjSplit{f}(u,d) = shapeLabel;
        end; clear d
    end; clear u
end; clear f


%% Make pie charts for direction and reward tuning
countPerDirRewTuning_bySubj = nan(nsubjects,4);
for f = 1:nsubjects
    sigConds = pval_Dir_Rew_subjSplit{f} < 0.05;
    tuningLabel = sigConds*[1 2 0]'; % 1 = dir only, 2 = rew only, 3 = both, 0 = none...
    tuningLabel(tuningLabel==0) = 4; % now 4 = none
    countPerDirRewTuning_bySubj(f,:) = histcounts(tuningLabel,0.5:4.5);%,'normalization','probability');
end; clear f

dirRewTuningColors = {...
    [0.9 0.1 0];   % direction only
    [0 0.1 0.9];   % reward only
    [0.7 0.1 0.7]; % both
    [0.8 0.8 0.8]; % none
    };
alpha = 0.7;

% flip for visualization
dirRewTuningColors = flip(dirRewTuningColors);
countPerDirRewTuning_bySubj = fliplr(countPerDirRewTuning_bySubj);

figure
for f = 1:nsubjects
    subplot(1,3,f); hold on
    p = pie(countPerDirRewTuning_bySubj(f,:));
    for t = 1:4
        p(t*2-1).FaceColor = dirRewTuningColors{t};
        p(t*2-1).FaceAlpha = alpha;
    end; clear t
    
    
    title(['Monkey ' subjectNames{f}(1)])
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    xticks([]); yticks([]);
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
end; clear f
set(gcf,'position',[4 737 769 211])

% Save it!
figname = ['FigX_singleUnitTuning_dirAndRew' ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Make pie charts splitting up reward tuning

% Get the fraction for each tuning per subject
ntunings = length(shapeNames);
countPerRewShape_bySubj = nan(nsubjects,ntunings);
for f = 1:nsubjects
    countPerRewShape_bySubj(f,:) = 100*histcounts(rewTuningShapeLabel_byDir_subjSplit{f}(:),...
        0.5:ntunings+0.5);%,'normalization','probability');
end; clear f

% Make the figure
alpha = 0.9;
rewTuningColors = {...
    [0 0.25 0.5],...
    [0 0.25 0.9],...
    [0 0.75 1],...
    [0 1 1],...
    [1 1 1],...
    [0.8 0.8 0.8]};

% flip for visualization
rewTuningColors = flip(rewTuningColors);
countPerRewShape_bySubj = fliplr(countPerRewShape_bySubj);

figure
for f = 1:nsubjects
    subplot(1,3,f); hold on
    p = pie(countPerRewShape_bySubj(f,:));
    for t = 1:ntunings
        p(t*2-1).FaceColor = rewTuningColors{t};
        p(t*2-1).FaceAlpha = alpha;
    end; clear t
    
    
    title(['Monkey ' subjectNames{f}(1)])
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    xticks([]); yticks([]);
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
end; clear f
set(gcf,'position',[4 440 769 211])

% Save it!
figname = ['FigX_singleUnitTuning_rewShape_byDir' ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Do this combined across dirs
alpha = 0.9;
rewTuningColors = {...
    [0.6 0 0],... % mono-up
    [0.85 0 0],... % S-L-up
    [0 0 0.6],... % mono-down
    [0 0 0.85],... % S-L down
    [0 0.7 0],... % U
    [0 0.9 0],... % inv-U
    [1 0.3 0],... % J-spec-up
    [0 0.3 1],... % J-spec-down
    [0.8 0.8 0.8]}; % none


% Get the fraction for each tuning per subject
shapeNames = {'mono-up','S-L-up','mono-down','S-L-down','U','inv-U','J-spec-up','J-spec-down','none'};
newOrder = [1 2 7 3 4 8 5 6 9]; % put the spec near the others of similar dir
ntunings = length(shapeNames);
countPerRewShape_bySubj = nan(nsubjects,ntunings);
for f = 1:nsubjects
    countPerRewShape_bySubj(f,:) = histcounts(rewTuningShapeLabel_subjSplit{f}(:),...
        0.5:ntunings+0.5);%,'normalization','probability');
end; clear f
shapeNames = shapeNames(newOrder);
countPerRewShape_bySubj = countPerRewShape_bySubj(:,newOrder);
rewTuningColors = rewTuningColors(newOrder);

% Make the figure
% flip for visualization
rewTuningColors = flip(rewTuningColors);
countPerRewShape_bySubj = fliplr(countPerRewShape_bySubj);

figure
for f = 1:nsubjects
    subplot(1,3,f); hold on
    p = pie(countPerRewShape_bySubj(f,:));
    for t = 1:ntunings
        p(t*2-1).FaceColor = rewTuningColors{t};
        p(t*2-1).FaceAlpha = alpha;
    end; clear t
    
    
%     title(['Monkey ' subjectNames{f}(1)])
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    xticks([]); yticks([]);
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
end; clear f
set(gcf,'position',[4 440 769 211])

% Save it!
figname = ['FigX_singleUnitTuning_rewShape' ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Get the total number of units tuned for each reward case
sigDirTuning = cellArrayToVector(cellfun(@(x) x(:,1)<0.05, pval_Dir_Rew_subjSplit,'uniformoutput',false));
dirCounts = histcounts(sigDirTuning);
disp([num2str(dirCounts(2)) '/' num2str(sum(dirCounts)) ' units tuned to direction (' ...
    num2str(round(dirCounts(2)/sum(dirCounts)*100,2)) '%)'])

sigRewTuning = cellArrayToVector(cellfun(@(x) x(:,2)<0.05, pval_Dir_Rew_subjSplit,'uniformoutput',false));
rewCounts = histcounts(sigRewTuning);
disp([num2str(rewCounts(2)) '/' num2str(sum(rewCounts)) ' units tuned to reward (' ...
    num2str(round(rewCounts(2)/sum(rewCounts)*100,2)) '%)'])

rewShapeCounts = histcounts(cellArrayToVector(rewTuningShapeLabel_byDir_subjSplit),0.5:ntunings+0.5);
tunedUnitDirCount = sum(rewShapeCounts(1:end-1));
disp([num2str(rewShapeCounts(1)) '/' num2str(sum(tunedUnitDirCount)) ' unit-dir show mono-up (' ...
    num2str(round(rewShapeCounts(1)/sum(tunedUnitDirCount)*100,2)) '%)'])
disp([num2str(rewShapeCounts(2)) '/' num2str(sum(tunedUnitDirCount)) ' unit-dir show mono-down (' ...
    num2str(round(rewShapeCounts(2)/sum(tunedUnitDirCount)*100,2)) '%)'])
disp([num2str(rewShapeCounts(3)) '/' num2str(sum(tunedUnitDirCount)) ' unit-dir show U (' ...
    num2str(round(rewShapeCounts(3)/sum(tunedUnitDirCount)*100,2)) '%)'])
disp([num2str(rewShapeCounts(4)) '/' num2str(sum(tunedUnitDirCount)) ' unit-dir show inv-U (' ...
    num2str(round(rewShapeCounts(4)/sum(tunedUnitDirCount)*100,2)) '%)'])


