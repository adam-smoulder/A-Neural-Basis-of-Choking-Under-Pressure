%%% This script makes Figure 1B, showing success rates as a function of
%%% reward.
%%%
%%% Adam Smoulder, 7/13/22

figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
filename_bySubject = {...
    'Fig1BehaviorData_Earl_20220816_134222';
    'Fig1BehaviorData_Prez_20220816_134302';
    'Fig1BehaviorData_Rocky_20220816_134525';
    };
nsubjects = length(filename_bySubject);
dateString = grabDateTimeString;

%% For each subject, get:
%  - overall mean success rate by reward
%  - error bars for overall by reward
%  - mean success rate for each day

subjectNames = cell(nsubjects,1);
successRate_subjSplit_byRew = cell(nsubjects,1);
pvals_subjSplit_byRewXRew = cell(nsubjects,1);
successRate_subjSplit_byRew_err = cell(nsubjects,1);
successRate_subjSplit_byRew_byDay = cell(nsubjects,1);
nboots = 10000; % number of bootstraps to use for errorbar SE
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
%     badInds = delayLengths < (400-8); % skip short delay inds
    if strcmp(taskInfo.subjectName,'Rocky') % only remove unprepared trials
        badInds = delayLengths < 317;
    else
        badInds = false(size(delayLengths));
    end
    
    % Group data by reward and get average/error
    ntrials = sum(~badInds);
    [~,~,n_byRew_bySF] = groupDataByLabel(ones(ntrials,1),... % 1st col = succ, 2nd = fail, rows = rewards
        [rewardLabels(~badInds) double(trialStatusLabels(~badInds)~=1)]);
    successRate_subjSplit_byRew{f} = 100*n_byRew_bySF(:,1)./sum(n_byRew_bySF,2);
    [~,~,~,~,~,~,SE] = bootstrapBinaryEvent(n_byRew_bySF(:,1),sum(n_byRew_bySF,2),nboots); % this is in counts, not %
    successRate_subjSplit_byRew_err{f} = 100*SE'./sum(n_byRew_bySF,2);
    
    % Get the p-value matrix from binomial proportion test
    [pvals_subjSplit_byRewXRew{f}] = binomialProportionTest2(n_byRew_bySF(:,1),sum(n_byRew_bySF,2));
    
    % Group data by reward x day and get average
    [~,~,n_byRew_byDay_bySF] = groupDataByLabel(ones(ntrials,1),... % 1st col = succ, 2nd = fail, rows = rewards
        [rewardLabels(~badInds) dayLabels(~badInds) double(trialStatusLabels(~badInds)~=1)]);
    successRate_subjSplit_byRew_byDay{f} = 100*n_byRew_byDay_bySF(:,:,1)./sum(n_byRew_byDay_bySF,3);
    
    sum(n_byRew_bySF,2)
end; clear f


%% Make a plot showing the individual days overlaid with the average
figure
indDayColor = [0.7 0.7 0.7];
meanColor = [0 0 0];
scatSize = 0.25; % how much scatter for each x to have for days
% minVal = floor(min(cellfun(@(x) min(x(:)),successRate_subjSplit_byRew_byDay))/5)*5-2; % round to nearest 5 then -1
% maxVal = round(max(cellfun(@(x) max(x(:)),successRate_subjSplit_byRew_byDay))/5)*5+2; % round to nearest 5 then +1
minVal = 42;
maxVal = 98;
rewColors = getDistinctColors('SELECT_ORDER',8);
rewNames = {'Small','Medium','Large','Jackpot'};
colorRewNames = rewNames;
for r = 1:nrewards
    colorRewNames{r} = ['\color{red} ' rewNames{r}];
    colorRewNames{r} = sprintf('\\color[rgb]{%f,%f,%f}%s', 0.9*rewColors{r}, rewNames{r});
end; clear r

for f = 1:nsubjects
    subplot(1,3,f); hold on
    ndays = size(successRate_subjSplit_byRew_byDay{f},2);
    xscats = linspace(-scatSize/2,scatSize/2,ndays);
    for a = 1:ndays
        plot(rewards+xscats(a), successRate_subjSplit_byRew_byDay{f}(:,a),'.-','markersize',5,'linewidth',0.5,'color',indDayColor)
    end; clear a
    errorbar(rewards,successRate_subjSplit_byRew{f},successRate_subjSplit_byRew_err{f},...
        '.-','markersize',2,'linewidth',2,'color',meanColor)
    
    axis([0.5 nrewards+0.5 minVal maxVal])
    xticks(1:4)
    xticklabels(colorRewNames)
    xtickangle(45)
    set(gca,'tickdir','out')
    title(['Monkey ' subjectNames{f}(1) ' (' num2str(ndays) ' sessions)'])
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    yticks([50 70 90])
    if f == 1
        ylabel('Success Rate (%)')
        yticklabels({'50' '70','90'})
    else
        yticklabels('')
    end
    
    % Show p-values for each S->L and L->J comparison
    disp(subjectNames{f})
    succRates = successRate_subjSplit_byRew{f};
    
    disp([num2str(diff(succRates([1 3]))) '%, p_{S->L} = ' num2str(pvals_subjSplit_byRewXRew{f}(1,3))])
    disp([num2str(diff(succRates([3 4]))) '%, p_{L->J} = ' num2str(pvals_subjSplit_byRewXRew{f}(3,4))])
end; clear n
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 724 432])

% Save the figure
figname = ['Fig1B_successRates']
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
