
filename_bySubject = {...
    ...'Fig1BehaviorData_Nelson_20221208_131310.mat';
    ...'Fig1BehaviorData_Ford_20221208_131312.mat';
    'Fig1BehaviorData_Earl_20220816_134222';
    'Fig1BehaviorData_Prez_20220816_134302';
    'Fig1BehaviorData_Rocky_20220816_134525';
    };
nsubjects = length(filename_bySubject);
figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;


%% For each animal, we'll get success rates as a function of delay length
nboots = 10000;
epochGroups = {[1 -11 -12 -22 -23 -34],[1 -22 -23 -34],[1 -34]}; % for consideration in delay, reach, and hold
epochSuccessGroups = {[1 -22 -23 -34],[1 -34],[1]};
delays_subjSplit = {[450:100:950],[250:200:850],[250:100:750]};
n_byRG_subjSplit = cell(nsubjects,1);
nfail_RG_subjSplit = cell(nsubjects,1);

for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = subjectName;
    delays = delays_subjSplit{f}; ndelays = length(delays);
    n_byRG_subjSplit{f} = zeros(4,ndelays);
    nfail_RG_subjSplit{f} = zeros(4,ndelays);
    for i = 1:length(trialStatusLabels)
        % oldchecksum = sum(n_byRG_subjSplit{f}(:));
        rewInd = rewardLabels(i);
        if ismember(trialStatusLabels(i),[1 -22 -23 -34]) % not delay failure
            [~, delayInd] = min(abs(delays - delayLengths(i)));
            n_byRG_subjSplit{f}(rewInd,delayInd) = 1 +  n_byRG_subjSplit{f}(rewInd,delayInd);
            if trialStatusLabels(i)~=1 % also add a fail
                nfail_RG_subjSplit{f}(rewInd,delayInd) = 1 +  nfail_RG_subjSplit{f}(rewInd,delayInd);
            end
        else % delay failure; find appropriate delays and distribute the failure
            validDelayInds = delays > (delayLengths(i)-17); % it could have been any of these delay lengths (give one sample back for cheats)
            if sum(validDelayInds) == 0 % likely a sample extra or Rocky long delay failure; in either case, only the last delay is valid
                validDelayInds(end) = true;
            end
            n_byRG_subjSplit{f}(rewInd,validDelayInds) = 1/sum(validDelayInds) +  n_byRG_subjSplit{f}(rewInd,validDelayInds);
            nfail_RG_subjSplit{f}(rewInd,validDelayInds) = 1/sum(validDelayInds) +  nfail_RG_subjSplit{f}(rewInd,validDelayInds);
        end
        % newchecksum = sum(n_byRG_subjSplit{f}(:));
        % assert(abs((newchecksum-oldchecksum)-1) < 1e-10)
    end; clear i
end; clear f

psucc_byRG_subjSplit = cellfun(@(x,y) (x-y)./x, n_byRG_subjSplit, nfail_RG_subjSplit, 'uniformoutput', false);
succRate_byRG_subjSplit = cellfun(@(x) 100*x, psucc_byRG_subjSplit, 'uniformoutput', false);
sem_byRG_subjSplit = cellfun(@(x,y) sqrt(x.*(1-x)./y)*100, psucc_byRG_subjSplit, n_byRG_subjSplit, 'uniformoutput', false);

%% Plot them
figure; set(gcf,'position',[-1511 304 1373 315])
rewColors = getDistinctColors('SELECT_ORDER',8);
for f = 1:nsubjects
    subplot(1,3,f); hold on
    x = delays_subjSplit{f};
    for r = 1:nrewards % For legend and lw purposes
        y = succRate_byRG_subjSplit{f}(r,:);
        plot(x,y,'-','linewidth',2,'color',rewColors{r})
    end; clear r
    for r = 1:nrewards
        y = succRate_byRG_subjSplit{f}(r,:);
        err = sem_byRG_subjSplit{f}(r,:);
        shadedErrorBar(x,y,err,'lineprops',{'color',rewColors{r}})
    end; clear r
    axis([-inf inf 30 95])
    if f==1
        ylabel('Success Rate (%)')
    end
    if f == 2
        xlabel('Delay duration (ms)')
    end
    title(['Monkey ' subjectNames{f}(1)])
    set(gca,'fontsize',12,'fontname','arial','tickdir','out')
end; clear f

% Save it!
figname = 'FigSX_successRateByDelayDuration'
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
