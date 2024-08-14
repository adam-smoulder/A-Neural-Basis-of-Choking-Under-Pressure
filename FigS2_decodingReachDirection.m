%%% This function tests how well we can decode the upcoming reach direction
%%% from a bin at the end of the delay period.
%%%
%%% This version trains decoders individually for each reward. To fairly
%%% compare results across rewards, we do normal cross validation with repeated
%%% subsamples. We do this because S/M/L rewards have many more trials than
%%% Jackpots, and hence, much more training data. We find the reward x
%%% direction condition with the least # of trials and subsample each
%%% condition down to that, then do K-fold cross-validation with the
%%% decoding, and store the average. We repeat this many times to
%%% approximate the average.
%%%
%%% Adam Smoulder, 7/21/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20240528_114423_newStitching';
    ...'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    ...'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;

rng('default'); rng(3195)

%% For each animal, get the decoding accuracy and error bars. Takes a while
%  to run.
methodNames = {'GNB','LDA','5NN'};
nmethods = length(methodNames);
withinReward = true; % true = fit one decoder for each reward, false = fit one decoder across rewards
nfolds = 5;
nreps_mean = 1000;
subjectNames = cell(nsubjects,1);
decodeAcc_byRew_byMethod_mean_subjSplit = cell(nsubjects,1);
n_byRew_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    
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
    
    % Do decoding using subsampled cross-validation
    decodeAcc_byRew_byMethod_mean = nan(nrewards,nmethods);
    for m = 1:nmethods
        disp(['Analyzing data for subject ' num2str(f) ', ' methodNames{m}])
        [decodeAcc_byRew_byMethod_mean(:,m), estSEM_byRew] = subsampledCV(...
            binnedData_pGC(~skipLabels,:),directionLabels(~skipLabels),rewardLabels(~skipLabels),...
            methodNames{m},nfolds,nreps_mean);
        disp('Completed subsampled CVing')
    end; clear m
    
    decodeAcc_byRew_byMethod_mean_subjSplit{f} = decodeAcc_byRew_byMethod_mean;
end; clear f

disp('Done processing data')



%% Plot results
axisLim_subjSplit = {[35 75], [70 90], [80 100]};
decoderColors = flip({[0 0.9 0],[0 0.55 0],[0 0.2 0]});

figure
rewColors = getDistinctColors('SELECT_ORDER',8);
rewNames = {'S','M','L','J'};
lw = 1.5;
ms = 15;
for f = 1:nsubjects
    subplot(1,3,f); hold on
    x = rewards;
    for m = 1:nmethods
        y = decodeAcc_byRew_byMethod_mean_subjSplit{f}(:,m);
        plot(x,y,'k-','linewidth',lw,'color',decoderColors{m})
    end; clear m
    for m = 1:nmethods
        y = decodeAcc_byRew_byMethod_mean_subjSplit{f}(:,m);
        for r = 1:nrewards
            plot(x(r), y(r),'k.-','linewidth',lw,'markersize',ms,'color',rewColors{r})
        end; clear r
    end; clear m
    
    % Labels and such
    set(gca,'fontsize',10,'fontname','arial','tickdir','out');
    xticks(rewards)
    xticklabels(rewNames)
    
    axis([min(rewards)-0.5 max(rewards)+0.5 axisLim_subjSplit{f}])
    title(['Monkey ' subjectNames{f}(1)])
    if f == 1
        ylabel('Direction Decode Accuracy (%)')
    elseif f == 2
        xlabel('Reward')
    end
end; clear f
subplot(1,nsubjects,2)
legend(methodNames,'location','s')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 495 300])

% % save it!
figname = ['FigS3_decodingReachDirection']
% saveFigAndSvg([figfolder 'recentFigs\'],figname);
% saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);




