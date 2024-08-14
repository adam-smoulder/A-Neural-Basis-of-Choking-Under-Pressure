%%% This script looks at single neuron tuning as a function of reward at
%%% the end of the delay period. We'll also look at directional tuning as
%%% well.
%%% 
%%% Adam Smoulder, 8/2/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
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
pval_DirTuning_subjSplit = cell(nsubjects,1);
rewTuningShapeLabel_subjSplit = cell(nsubjects,1);
pvalDir_byRew_subjSplit = cell(nsubjects,1);
rewTuningShapeLabel_byDir_subjSplit = cell(nsubjects,1);
rewTuningChanges_subjSplit = cell(nsubjects,1);
projDataMeans_byDir_byRew_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    [ntrials,~] = size(singleNeuronData_pGC);
    
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
    
    singleNeuronData_pGC(skipLabels,:) = [];
    rewardLabels(skipLabels) = [];
    directionLabels(skipLabels) = [];
    
    % Identify units with at least nmin of each direction x reward
    % condition and throw away all else. 
    [validUnitInds,nminByUnit] = unitsThatHaveMinTrialCount(singleNeuronData_pGC,...
        nmin_bySubj(f),[directionLabels rewardLabels]);
    nunits = sum(validUnitInds);
    validNeuronData = singleNeuronData_pGC(:,validUnitInds);

    % For the remaining valid units, get tuning properties
    for u = 1:nunits
        curData = validNeuronData(:,u);
        validInds = ~isnan(curData);
        meanUnitFR_byDir_byRew_subjSplit{f}(u,:,:) = cellfun(@(x) mean(x), groupDataByLabel(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)]));
        semUnitFR_byDir_byRew_subjSplit{f}(u,:,:) = cellfun(@(x) nansem(x,1), groupDataByLabel(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)]));
        varUnitFR_byDir_byRew_subjSplit{f}(u,:,:) = cellfun(@(x) var(x), groupDataByLabel(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)]));
    
        [pvals,~,stats] = anovan(curData(validInds),[directionLabels(validInds) rewardLabels(validInds)],'display','off','model','interaction');
        pval_DirTuning_subjSplit{f}(u) = pvals(1);
    end; clear u
    
    % Do PCA on the reward averages to get reward axis
    rewMeans = squeeze(mean(meanUnitFR_byDir_byRew_subjSplit{f},2))';
    mu = mean(rewMeans);
    [wR,zR,~,~,pveR] = pca(rewMeans-mu,'numcomponents',1);
    if zR(1) > zR(3) % flip the sign so L > S
        wR = -wR; zR = -zR;
    end
        
    % Do PCA on the target averages to get target plane
    dirMeans = mean(meanUnitFR_byDir_byRew_subjSplit{f},3)';
    [wD,zD,~,~,pveD] = pca(dirMeans-mu,'numcomponents',2);
    
    subspace(wR,wD)*180/pi
    
    % Project all data down
    [w,~] = qr([wD wR]);
    w = w(:,1:3);
    if rewMeans(1,:)*w(:,3) > rewMeans(3,:)*w(:,3)
        w(:,3) = -w(:,3);
    end
    projDataMeans_byDir_byRew = cell(ndirections,nrewards);
    for d = 1:ndirections
        for r = 1:nrewards
            projDataMeans_byDir_byRew{d,r} = (meanUnitFR_byDir_byRew_subjSplit{f}(:,d,r)'-mu)*w;
        end; clear r
    end; clear d
    projDataMeans_byDir_byRew_subjSplit{f} = projDataMeans_byDir_byRew;
    
    % Get fano factor
    FFUnitFR_byDir_byRew_subjSplit{f} = varUnitFR_byDir_byRew_subjSplit{f} ./ (meanUnitFR_byDir_byRew_subjSplit{f}+1);
end; clear f


%% Show a single example TC
f = 1;
u = 27;
r = 2;
y = squeeze(meanUnitFR_byDir_byRew_subjSplit{f}(u,:,r));

b0 = squeeze(mean(y,2));
b1 = squeeze(2/8*((y(:,2,:)+y(:,4,:)-y(:,6,:)-y(:,8,:))/sqrt(2)+(y(:,3,:)-y(:,7,:))));
b2 = squeeze(2/8*((y(:,2,:)-y(:,4,:)-y(:,6,:)+y(:,8,:))/sqrt(2)+(y(:,1,:)-y(:,5,:))));

x = 0:45:315;
x_cont = linspace(-22.5,360-22.5,100);
y_cont = b0 + b1*sind(x_cont) + b2*cosd(x_cont);

figure; hold on
lw = 2.5;
ms = 45;
rewColors = getDistinctColors('SELECT_ORDER',8);
plot(x_cont,y_cont,'-','color',rewColors{r},'linewidth',lw)
plot(x,y,'.','color',0.9*rewColors{r},'markersize',ms)
plot(x_cont,b0*ones(length(x_cont),1),'k-','linewidth',0.5)
axis([min(x_cont) max(x_cont) 19 41])
set(gca,'fontsize',14,'fontname','arial','tickdir','out','box','off');
xlabel(['Reach direction'])
xticks(x)
xticklabels({['0' char(176)],['45' char(176)],['90' char(176)],['135' char(176)],['180' char(176)],['225' char(176)],['270' char(176)],['315' char(176)]})
ylabel('Firing Rate (sp/s)')

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 547 314])

% Save it!
figname = ['FigS4_singleUnitModDepth_ExampleTC_r' num2str(r) '_f' num2str(f) 'u' num2str(u)]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Alternatively, show each reward for one neuron
f = 1;
u = 36;
meanFRs = squeeze(meanUnitFR_byDir_byRew_subjSplit{f}(u,:,:));

figure
lw = 2.5;
ms = 35;
rewColors = getDistinctColors('SELECT_ORDER',8);
for r = 1:nrewards
    y = meanFRs(:,r)';
    b0 = squeeze(mean(y,2));
    b1 = squeeze(2/8*((y(:,2,:)+y(:,4,:)-y(:,6,:)-y(:,8,:))/sqrt(2)+(y(:,3,:)-y(:,7,:))));
    b2 = squeeze(2/8*((y(:,2,:)-y(:,4,:)-y(:,6,:)+y(:,8,:))/sqrt(2)+(y(:,1,:)-y(:,5,:))));
    
    x = 0:45:315;
    x_cont = linspace(-22.5,360-22.5,100);
    y_cont = b0 + b1*sind(x_cont) + b2*cosd(x_cont);
    
    subplot(1,nrewards,r); hold on
    plot(x_cont,y_cont,'-','color',rewColors{r},'linewidth',lw)
    plot(x,y,'.','color',0.9*rewColors{r},'markersize',ms)
    plot(x_cont,b0*ones(length(x_cont),1),'-','linewidth',0.5,'color',rewColors{r})
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    xticks(x)
    xticklabels({['0' char(176)],['45' char(176)],['90' char(176)],['135' char(176)],['180' char(176)],['225' char(176)],['270' char(176)],['315' char(176)]})
    if u == 27
        yticks([20 30 40])
        axis([min(x_cont) max(x_cont) 19 41])
    elseif u == 36
        yticks([0 10 20])
        axis([min(x_cont) max(x_cont) 0 22])
    else
       
        axis([min(x_cont) max(x_cont) min(meanFRs(:)) max(meanFRs(:))])
    end
    if r == 1
        xlabel(['Reach direction'])
        ylabel('Firing Rate (sp/s)')
    else
        yticklabels([])
    end
    
end; clear r

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 1665 314])

% Save it!
figname = ['FigS4_singleUnitModDepth_ExampleTCs_allRew_f' num2str(f) 'u' num2str(u) ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

%% Alternatively, show each reward for one neuron
f = 1;
u = 36;
meanFRs = squeeze(meanUnitFR_byDir_byRew_subjSplit{f}(u,:,:));
semFRs = squeeze(semUnitFR_byDir_byRew_subjSplit{f}(u,:,:));


figure
lw = 2.5;
ms = 25;
rewColors = getDistinctColors('SELECT_ORDER',8);
for r = 1:nrewards
    y = meanFRs(:,r)';
    b0 = squeeze(mean(y,2));
    b1 = squeeze(2/8*((y(:,2,:)+y(:,4,:)-y(:,6,:)-y(:,8,:))/sqrt(2)+(y(:,3,:)-y(:,7,:))));
    b2 = squeeze(2/8*((y(:,2,:)-y(:,4,:)-y(:,6,:)+y(:,8,:))/sqrt(2)+(y(:,1,:)-y(:,5,:))));
    
    x = 0:45:315;
    x_cont = linspace(-22.5,360-22.5,100);
    y_cont = b0 + b1*sind(x_cont) + b2*cosd(x_cont);
    
    err = semFRs(:,r)';
    
    subplot(1,nrewards,r); hold on
    plot(x_cont,y_cont,'-','color',rewColors{r},'linewidth',lw)
%     plot(x,y,'.','color',0.9*rewColors{r},'markersize',ms)
    errorbar(x,y,err,'.','color',0.9*rewColors{r},'markersize',ms)
    plot(x_cont,b0*ones(length(x_cont),1),'-','linewidth',0.5,'color',rewColors{r})
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    xticks(x)
    xticklabels({['0' char(176)],['45' char(176)],['90' char(176)],['135' char(176)],['180' char(176)],['225' char(176)],['270' char(176)],['315' char(176)]})
    if u == 27
        yticks([20 30 40])
        axis([min(x_cont) max(x_cont) 19 41])
    elseif u == 36
        yticks([0 10 20])
        axis([min(x_cont) max(x_cont) 0 25])
    else
       
        axis([min(x_cont) max(x_cont) min(meanFRs(:)) max(meanFRs(:))])
    end
    if r == 1
        xlabel(['Reach direction'])
        ylabel('Firing Rate (sp/s)')
    else
        yticklabels([])
    end
    
end; clear r

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 1665 314])

% Save it!
figname = ['FigS4_singleUnitModDepth_ExampleTCsErrBar_allRew_f' num2str(f) 'u' num2str(u) ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

%% Show modulation depth as a function of reward
nunitsToPlot = 20;
lw = 2;
ms = 13;
figure
rewColors = getDistinctColors('SELECT_ORDER',8);
meanUnitFR_byDR_sorted_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    subplot(nsubjects,5,(1:4)+(f-1)*5)
    nunits = length(pvalDir_byRew_subjSplit{f});
    y = meanUnitFR_byDir_byRew_subjSplit{f};
    
    % From Georgopoulos et al 1982, we get coefficients for the form of
    % y = b0 + b1*sin(theta) + b2*cos(theta)
    if size(y,2)==8 % 8 directions
        b0 = squeeze(mean(y,2));
        b1 = squeeze(2/8*((y(:,2,:)+y(:,4,:)-y(:,6,:)-y(:,8,:))/sqrt(2)+(y(:,3,:)-y(:,7,:))));
        b2 = squeeze(2/8*((y(:,2,:)-y(:,4,:)-y(:,6,:)+y(:,8,:))/sqrt(2)+(y(:,1,:)-y(:,5,:))));
    elseif size(y,2)==4 % 4 directions
        b0 = squeeze(mean(y,2));
        b1 = squeeze((1/(sqrt(2)*2))*(y(:,1,:)+y(:,2,:)-y(:,3,:)-y(:,4,:)));
        b2 = squeeze((1/(sqrt(2)*2))*(y(:,1,:)-y(:,2,:)-y(:,3,:)+y(:,4,:)));
    end
    % This is actually pretty easy to derive by hand too!
    
    % Then, modulation depth (m) is the norm of these
    m = sqrt(b1.^2+b2.^2);
    [~,order] = sort(m(:,2),'descend');
    m = m(order,:);
    m_bySubject{f} = m;
    pd = mod(atan2d(b2,b1)+360,360);
    pd_bySubject{f} = pd(order,:);
    b0_bySubject{f} = b0(order,:);
    meanUnitFR_byDR_sorted_subjSplit{f} = meanUnitFR_byDir_byRew_subjSplit{f}(order,:,:);
    for u = 1:nunitsToPlot
        x = (1:nrewards)+(u-1)*(nrewards+1);
        y = m(u,:);
        plot(x,y,'k-','color',[0 0.9 0.9],'linewidth',lw); hold on
        %         semilogy(x,y,'k-','color',[0 0.9 0.9],'linewidth',lw); hold on
        for r = 1:nrewards
            plot(x(r),y(r),'.','markersize',ms,'color',rewColors{r})
            %             semilogy(x(r),y(r),'.','markersize',ms,'color',rewColors{r})
        end; clear r
    end; clear u
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    if f == 3
        xlabel('Unit # (sorted by mod. depth)')
    end
    
    ylabel(['Monkey ' subjectNames{f}(1)])
    if f == 1
        title(['Modulation depth of top ' num2str(nunitsToPlot) ' units'])
    end
%     xticks([])
    xticks(2.5:5:100)
    if f ==3
        xticklabels(((2.5:5:100)+2.5)/5)
    else
        xticklabels([])
    end
    
    % Show the average difference in modulation depth for the top 1/4 of units
    subplot(nsubjects,5,f*5); hold on
    fracNum = 4;
    n = floor(length(m)/fracNum)
%     n = 20; % Steve recommended we just use the units shown... but this doesn't make any sense. That's just sensitive to sampling error
    y = m(1:n,:)-mean(m(1:n,:),2); % de-mean for comparison
%     y = m(end-n:end,:)-mean(m(end-n:end,:),2); % de-mean for comparison
%     x = 0.1*randn(n,4)+(1:4);
%     for u = 1:length(x)
%         plot(x(u,:),y(u,:),'k-','linewidth',0.5)
%     end; clear u
    y_mean = mean(y)+mean(mean(m(1:n,:),2)); % Add overall mean of mod-depth for these for reference
    err = nansem(y,1);
    minVal = min(y_mean(:)-err(:))-0.1*range(y_mean(:));
    maxVal = max(y_mean(:)+err(:))+0.1*range(y_mean(:));
%     errorbar(y_mean,err,'k-','linewidth',lw,'color',[0 0.9 0.9])
    plot(y_mean,'k-','linewidth',lw,'color',[0 0.9 0.9])
    for r = 1:nrewards
        errorbar(r,y_mean(r),err(r),'k.-','linewidth',lw,'color',rewColors{r},'markersize',2*ms)
    end; clear r
    axis([0.5 4.5 minVal maxVal])
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    if f == 1
        title(['Avg. of top 1/' num2str(fracNum) ' of units'])
    end
    xticks(1:4)
    if f == 3
        xticklabels({'S','M','L','J'})
    else
        xticklabels([])
    end
    
    % Try stats
    [~,pSL(f)] = ttest2(y(:,1),y(:,3));
    [~,pLJ(f)] = ttest2(y(:,3),y(:,4));
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[0 500] 1369 933])
pSL
pLJ

% Save the figure
% Save it!
figname = ['FigS4_singleUnitModDepth_fracNum' num2str(fracNum) ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Plot the rings
figure
set(gcf,'renderer','painters')
rewColors = getDistinctColors('SELECT_ORDER',8);
dirColors = getDistinctColors('SELECT_ORDER',16);
ms1 = 9;
lw1 = 3;
lw2 = 3;

view_subjSplit = {[-30 15],[-140 10],[80 20]}; % good for PCA axis equal
tickVals_subjSplit = {[-210:35:210],[-210:35:210],[-240:60:240]};
axis_subjSplit = {[-70 105 -105 70 -70 70],[-75 75 -105 105 -70 105],[-180 180 -240 180 -75 50]};


for f = 1:nsubjects
    % Get the current points
    if strcmp(subjectNames{f},'Prez')
        directions = 2:2:8;
    else
        directions = 1:8;
    end
    ndirections = length(directions);
    
    % Show the reward axis as well
    subplot(1,nsubjects,f); hold on
    rewTargProj_byDir_byRew_means = projDataMeans_byDir_byRew_subjSplit{f};
    for r = 1:nrewards
        curMeans = reshape([rewTargProj_byDir_byRew_means{:,r}],[3 ndirections])';
        curMeans(:,1:2) = curMeans(:,1:2)-mean(curMeans(:,1:2)); % we can center within rewards when projecting in target plane
        ringPoints = [curMeans ; curMeans(1,:)];
        plot3(ringPoints(:,1),ringPoints(:,2),ringPoints(:,3),'-','linewidth',1.5*lw1,...
            'color',rewColors{r})
        for d = 1:ndirections
            plot3(curMeans(d,1),curMeans(d,2),curMeans(d,3),'o','markersize',1.25*ms1,'linewidth',lw1,...
                'markeredgecolor',rewColors{rewards(r)},'markerfacecolor',dirColors{directions(d)});
        end; clear d
    end; clear r
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    if f == 1
        zlabel('Unit Reward Axis')
    end
    %     xlabel('Unit Targ. Ax. 1')
    %     ylabel('Unit Targ. Ax. 2')
    xlabel('Ax. 1') % makes it look better; fix in illustrator
    ylabel('Ax. 2')
    title(['Monkey ' subjectNames{f}(1)])
    view(view_subjSplit{f})
    axis equal
    grid on
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[0 500] 1680 933])

% Save it
figname = ['FigS4_singleUnitRewardRings']
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Show fano factor for each reach direction 
figure; set(gcf,'position',[2123 565 844 353])
for f = 1:nsubjects
    % Concatenate across dirs
    curData = permute(FFUnitFR_byDir_byRew_subjSplit{f},[3 1 2]);
    curData = curData(:,:)';
%     curData = curData - mean(curData,2); % subtract avg over rewards for each unit x dir
    
    
    subplot(1,nsubjects,f); hold on
    x = 1:4;
    y = mean(curData);
    err = nansem(curData,1);
    errorbar(x,y,err)
    axis([0.5 4.5 -inf inf])
    title({['Monkey ' subjectNames{f}(1)],['n = ' num2str(size(curData,1)) ' units x directions']})
    set(gca,'fontsize',10,'fontname','arial','tickdir','out');
    ylabel('Fano Factor (sp/s)')
end; clear f

figure; set(gcf,'position',[2123 565 844 353])
for f = 1:nsubjects
    % Concatenate across dirs
    curData = permute(FFUnitFR_byDir_byRew_subjSplit{f},[3 1 2]);
    curData = curData(:,:)';
    curData = curData - mean(curData,2); % subtract avg over rewards for each unit x dir
    
    
    subplot(1,nsubjects,f); hold on
    x = 1:4;
    y = mean(curData);
    err = nansem(curData,1);
    errorbar(x,y,err)
    axis([0.5 4.5 -inf inf])
    title({['Monkey ' subjectNames{f}(1)],['n = ' num2str(size(curData,1)) ' units x directions']})
    set(gca,'fontsize',10,'fontname','arial','tickdir','out');
    ylabel('Fano Factor (unit-dir avg = 0) (sp/s)')
end; clear f

%% Show noise variance as a function of reward
nunitsToPlot = 20;
lw = 2;
ms = 10;
figure
rewColors = getDistinctColors('SELECT_ORDER',8);
meanUnitFR_byDR_sorted_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    subplot(nsubjects,5,(1:4)+(f-1)*5)
    nunits = length(pvalDir_byRew_subjSplit{f});
    noiseVar = squeeze(mean(varUnitFR_byDir_byRew_subjSplit{f},2));
    sigVar = squeeze(var(meanUnitFR_byDir_byRew_subjSplit{f},[],2));

    % Sort by (biased) signal variance
    [~,order] = sort(sigVar(:,2),'descend');

%     % Sort by noise
%     [~,order] = sort(noiseVar(:,2),'descend');

%     % Sort by SNR
%     [~,order] = sort(sigVar(:,2)./noiseVar(:,2),'descend');

    noiseVar = noiseVar(order,:);
    noiseVar_bySubject{f} = m;
    for u = 1:nunitsToPlot
        x = (1:nrewards)+(u-1)*(nrewards+1);
        y = sqrt(noiseVar(u,:));
        plot(x,y,'k-','color',[0.9 0 0.9],'linewidth',lw); hold on
        %         semilogy(x,y,'k-','color',[0 0.9 0.9],'linewidth',lw); hold on
        for r = 1:nrewards
            plot(x(r),y(r),'.','markersize',ms,'color',rewColors{r})
            %             semilogy(x(r),y(r),'.','markersize',ms,'color',rewColors{r})
        end; clear r
    end; clear u
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    if f == 3
%         xlabel('Unit # (sorted by mod. depth)')
%         xlabel('Unit # (sorted by noise var)')
        xlabel('Unit # (sorted by SNR)')
    end

    ylabel(['Monkey ' subjectNames{f}(1)])
    if f == 1
        title(['Noise (Hz) of top ' num2str(nunitsToPlot) ' units'])
    end
%     xticks([])
    xticks(2.5:5:100)
    if f ==3
        xticklabels(((2.5:5:100)+2.5)/5)
    else
        xticklabels([])
    end

    % Show the average difference in modulation depth for the top 1/4 of units
    subplot(nsubjects,5,f*5); hold on
    fracNum = 1;
    n = floor(length(noiseVar)/fracNum);
    y = noiseVar(1:n,:)-mean(noiseVar(1:n,:),2); % de-mean for comparison
    y_mean = mean(y);
    err = nansem(y,1);
    errorbar(y_mean,err,'k-','linewidth',lw)
    axis([0.5 4.5 -inf inf])
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    if f == 1
        title(['Avg. of top 1/' num2str(fracNum) ' of units'])
    end
    xticks(1:4)
    if f == 3
        xticklabels({'S','M','L','J'})
    else
        xticklabels([])
    end

    % Try stats
    [~,pSL(f)] = ttest2(y(:,1),y(:,3));
    [~,pLJ(f)] = ttest2(y(:,3),y(:,4));
end; clear f
set(gcf,'position',[1 31 1369 933])
pSL
pLJ

% Save the figure
% Save it!
figname = ['FigX_singleUnitTuning_noiseVar_fracNum' num2str(fracNum) ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

% 
% %% Look at avg firing rate versus modulation depth
% figure
% set(gcf,'position',[13 528 1661 420])
% ms = 10;
% for f = 1:nsubjects
%     subplot(1,3,f); hold on
%     x = b0_bySubject{f};
%     y = m_bySubject{f};
%     
%     % Show only best tuned units
%     nToPlot = 5; %floor(length(y)/4);
%     [~,order] = sort(y(:,2),'descend');
%     x = x(order(1:nToPlot),:);
%     y = y(order(1:nToPlot),:);
%     plot(x',y','k-')
%     for r = 1:nrewards
%         plot(x(:,r),y(:,r),'.','markersize',ms,'color',rewColors{r})
%     end; clear r
%     corr(x(:),y(:),'type','spearman')
%     
% %     plot(mean(x), mean(y),'linewidth',5)
%     
%     xlabel('Mean FR (sp/s)')
%     ylabel('Mod. Depth (sp/s)')
%     set(gca,'fontsize',14,'fontname','arial','tickdir','out','box','off');
% %     axis equal
% end; clear f
% 
% 
% %% Show the same idea, but with difference from Small to Large
% figure
% set(gcf,'position',[13 528 1661 420])
% ms = 10;
% for f = 1:nsubjects
%     subplot(1,3,f); hold on
%     pvals = pval_DirTuning_subjSplit{f};
%     x = diff(b0_bySubject{f}(pvals<0.05/length(x),[1 3]),[],2);
%     y = diff(m_bySubject{f}(pvals<0.05/length(x),[1 3]),[],2);
% %     size(x) 
%     
%     % Only do this with the best tuned units
%     nToShow = length(x)/4;%length(y)/4;%20; %floor(length(y)/4);
%     [~,order] = sort(m_bySubject{f}(:,2),'descend');
%     x = x(order(1:nToShow));
%     y = y(order(1:nToShow));
%     %
%     plot(linspace(-max(abs(x)),max(abs(x)),100),zeros(100,1),'-','color',0.8*ones(3,1))
%     plot(zeros(100,1),linspace(-max(abs(y)),max(abs(y)),100),'-','color',0.8*ones(3,1))
%     plot(x,y,'.','markersize',ms,'color',rewColors{1})
%     plot(mean(x),mean(y),'^','markersize',ms,'linewidth',2,'color',0.8*rewColors{1})
%     corr(x(:),y(:),'type','spearman')
%     axis([-max(abs(x)) max(abs(x)) -max(abs(y)) max(abs(y))])
% 
% %     plot(mean(x), mean(y),'linewidth',5)
%     title(['Monkey ' subjectNames{f}(1)])
%     xlabel('\Delta_{S,L} Mean FR (sp/s)')
%     ylabel('\Delta_{S,L} Mod. Depth (sp/s)')
%     set(gca,'fontsize',14,'fontname','arial','tickdir','out','box','off');
% %     axis equal
% end; clear f
% 
% 
% %% Show the same idea, but with difference from Large to Jackpot
% figure
% set(gcf,'position',[13 528 1661 420])
% ms = 10;
% for f = 1:nsubjects
%     subplot(1,3,f); hold on
%     pvals = pval_DirTuning_subjSplit{f};
%     x = diff(b0_bySubject{f}(pvals<0.05/length(x),[3 4]),[],2);
%     y = diff(m_bySubject{f}(pvals<0.05/length(x),[3 4]),[],2);
%     
%     % Only do this with the best tuned units
%     nToShow = length(x)/4; %floor(length(y)/4);
%     [~,order] = sort(m_bySubject{f}(:,2),'descend');
%     x = x(order(1:nToShow));
%     y = y(order(1:nToShow));
%     %
%     plot(linspace(-max(abs(x)),max(abs(x)),100),zeros(100,1),'-','color',0.8*ones(3,1))
%     plot(zeros(100,1),linspace(-max(abs(y)),max(abs(y)),100),'-','color',0.8*ones(3,1))
%     plot(x,y,'.','markersize',ms,'color',rewColors{4})
%     plot(mean(x),mean(y),'^','markersize',ms,'linewidth',2,'color',0.8*rewColors{4})
%     corr(x(:),y(:),'type','spearman')
%     axis([-max(abs(x)) max(abs(x)) -max(abs(y)) max(abs(y))])
% 
% %     plot(mean(x), mean(y),'linewidth',5)
%     title(['Monkey ' subjectNames{f}(1)])
%     xlabel('\Delta_{L,J} Mean FR (sp/s)')
%     ylabel('\Delta_{L,J} Mod. Depth (sp/s)')
%     set(gca,'fontsize',14,'fontname','arial','tickdir','out','box','off');
% %     axis equal
% end; clear f
% 
% %% Same as last time, but flip sign if S->L is negative
% figure
% set(gcf,'position',[13 528 1661 420])
% ms = 10;
% for f = 1:nsubjects
%     subplot(1,3,f); hold on
%     pvals = pval_DirTuning_subjSplit{f};
%     x = diff(b0_bySubject{f}(pvals < 1E-10/length(pvals),[3 4]),[],2);
%     y = diff(m_bySubject{f}(pvals < 1E-10/length(pvals),[3 4]),[],2);
%     
%     % If S->L is negative, flip sign for x
%     x = x.*sign(diff(b0_bySubject{f}(pvals < 1E-10/length(pvals),[1 3]),[],2));
%     
%     % Only do this with the best tuned units
%     nToShow = length(y); %floor(length(y)/4);
%     [~,order] = sort(m_bySubject{f}(:,2),'descend');
%     x = x(order(1:nToShow));
%     y = y(order(1:nToShow));
%     %
%     plot(linspace(-max(abs(x)),max(abs(x)),100),zeros(100,1),'-','color',0.8*ones(3,1))
%     plot(zeros(100,1),linspace(-max(abs(y)),max(abs(y)),100),'-','color',0.8*ones(3,1))
%     plot(x,y,'.','markersize',ms,'color',rewColors{4})
%     plot(mean(x),mean(y),'^','markersize',ms,'linewidth',2,'color',0.8*rewColors{4})
%     corr(x(:),y(:),'type','spearman')
%     axis([-max(abs(x)) max(abs(x)) -max(abs(y)) max(abs(y))])
% 
% %     plot(mean(x), mean(y),'linewidth',5)
%     title(['Monkey ' subjectNames{f}(1)])
%     xlabel('\Delta_{L,J} Mean FR (sp/s)')
%     ylabel('\Delta_{L,J} Mod. Depth (sp/s)')
%     set(gca,'fontsize',14,'fontname','arial','tickdir','out','box','off');
% %     axis equal
% end; clear f


% %% Look at diff in PD from med reward
% f = 3;
% curPDs = pd_bySubject{f};
% medDiffCenterPDs = (curPDs-curPDs(:,4))';
% % medDiffCenterPDs = (curPDs-mean(curPDs,2))';
% temp = min(abs(medDiffCenterPDs),abs(360+medDiffCenterPDs));
% minDiffCenterPDErrs = min(temp,360-temp);
% 
% nanmedian(minDiffCenterPDErrs')
% nansemed(minDiffCenterPDErrs',1)
% 
% figure; hold on
% plot(medDiffCenterPDs,'.-','color',0.8*[1 1 1],'linewidth',0.5)
% plot(mean(medDiffCenterPDs'),'k.-','linewidth',2)


% %% Now to test clipping, look at highest vs. lowest FR for neurons
% 
% 
% % ADAM HERE maybe look at TCs a little more to evaluate tuning shapes...?


%%


%% Show fano factor as a function of reward
nunitsToPlot = 20;
lw = 2;
ms = 13;
figure
rewColors = getDistinctColors('SELECT_ORDER',8);
meanUnitFR_byDR_sorted_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    subplot(nsubjects,5,(1:4)+(f-1)*5)
    nunits = length(pvalDir_byRew_subjSplit{f});
    sigVar = squeeze(var(meanUnitFR_byDir_byRew_subjSplit{f},[],2));
    FFs = squeeze(mean(FFUnitFR_byDir_byRew_subjSplit{f},2));

    % Sort by (biased) signal variance
    [~,order] = sort(sigVar(:,2),'descend');

%     % Sort by noise
%     [~,order] = sort(noiseVar(:,2),'descend');

%     % Sort by SNR
%     [~,order] = sort(sigVar(:,2)./noiseVar(:,2),'descend');

    FFs = FFs(order,:);
    for u = 1:nunitsToPlot
        x = (1:nrewards)+(u-1)*(nrewards+1);
        y = sqrt(FFs(u,:));
        plot(x,y,'k-','color',[0.9 0 0.9],'linewidth',lw); hold on
        %         semilogy(x,y,'k-','color',[0 0.9 0.9],'linewidth',lw); hold on
        for r = 1:nrewards
            plot(x(r),y(r),'.','markersize',ms,'color',rewColors{r})
            %             semilogy(x(r),y(r),'.','markersize',ms,'color',rewColors{r})
        end; clear r
    end; clear u
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    if f == 3
%         xlabel('Unit # (sorted by mod. depth)')
%         xlabel('Unit # (sorted by noise var)')
        xlabel('Unit # (sorted by mod. depth)')
    end

    ylabel(['Monkey ' subjectNames{f}(1)])
    if f == 1
        title(['Fano Factor (sp / s) of top ' num2str(nunitsToPlot) ' units'])
    end
%     xticks([])
    xticks(2.5:5:100)
    if f ==3
        xticklabels(((2.5:5:100)+2.5)/5)
    else
        xticklabels([])
    end

    % Show the average difference in modulation depth for the top 1/4 of units
    subplot(nsubjects,5,f*5); hold on
    fracNum = 1;
    n = floor(length(FFs)/fracNum);
    unitMeans = mean(FFs(1:n,:),2);
    grandMean = mean(unitMeans); % used to re-center
    y = FFs(1:n,:) - unitMeans + grandMean; % de-mean for comparison
    y_mean = mean(y);
    err = nansem(y,1);
    plot(y_mean,'k-','linewidth',lw,'color',[0.9 0 0.9])
    for r = 1:nrewards
        errorbar(r,y_mean(r),err(r),'k.-','linewidth',lw,'color',rewColors{r},'markersize',2*ms)
    end; clear r
    axis([0.5 4.5 -inf inf])
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    if f == 1
        % title(['Avg. of top 1/' num2str(fracNum) ' of units'])
        title('Fano Factor (sp/s)')
    end
    xticks(1:4)
    if f == 3
        xticklabels({'S','M','L','J'})
    else
        xticklabels([])
    end

    % Try stats
    [~,pSL(f)] = ttest2(y(:,1),y(:,3));
    [~,pLJ(f)] = ttest2(y(:,3),y(:,4));
end; clear f
set(gcf,'position',[1 31 1369 933])
pSL
pLJ


% Save it!
figname = ['FigS4_singleUnitFanoFactor_fracNum' num2str(fracNum) ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
