% Load PSTH data and make plots for it
% Adam Smoulder, 10/26/23

figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;
load('Fig1PSTHData_Earl_20231026_092738')
nrewards = 4;
binSize = 5;

%% Plot it!
unitsToPlot = [133 34 119];

figure
rewColors = getDistinctColors('SELECT_ORDER',8);
lw = 2;
indsToPlot_delay = find(time_delay >= -170 & time_delay <=500);
indsToPlot_GC = time_GC >= -200 & time_GC <= 100;
for u = 1:length(unitsToPlot)
    subplot(length(unitsToPlot),3,3*(u-1)+(1:2)); hold on
    unit = unitsToPlot(u);
    maxVal = -inf; minVal = inf;
    for r = 1:nrewards
        x = time_delay(indsToPlot_delay);
        y = delayTraj_byRew_nanmean{r}(unit,indsToPlot_delay);
        err = delayTraj_byRew_nansem{r}(unit,indsToPlot_delay);
        shadedErrorBar(x,y,err,...
            'lineProps',{'r-','color',rewColors{r},'linewidth',lw})
        maxVal = max(y(:)+err(:),maxVal);
        minVal = min(y(:)-err(:),minVal);
    end; clear r
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    if u ~= length(unitsToPlot)
        xticklabels([])
    else
        xlabel('Time from Target Onset (ms)')
    end
    ylabel(['Unit ' num2str(unit) ' FR'])
    axis([-inf inf -inf inf])
    
    subplot(length(unitsToPlot),3,3*u); hold on
    for r = 1:nrewards
        x = time_GC(indsToPlot_GC);
        y = GCTraj_byRew_nanmean{r}(unit,indsToPlot_GC);
        err = GCTraj_byRew_nansem{r}(unit,indsToPlot_GC);
        shadedErrorBar(x,y,err,...
            'lineProps',{'r-','color',rewColors{r},'linewidth',lw})
        maxVal = max([y(:)+err(:) ; maxVal]);
        minVal = min([y(:)-err(:) ; minVal]);
    end; clear r
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    if u ~= length(unitsToPlot)
        xticklabels([])
    else
        xlabel('Time from GC (ms)')
    end
    yticklabels([])
    
    % Round axes to nearest 5 and set them
    maxVal = ceil(maxVal/5)*5+1;
    minVal = floor(minVal/5)*5-1;
    axis([-inf inf minVal maxVal])
    
    subplot(length(unitsToPlot),3,3*(u-1)+(1:2)); hold on
    axis([-inf inf minVal maxVal])

end; clear u
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 716 443])

% Save it!
figname = ['Fig1C_PSTHExamples_labeled']
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Plot it but with much less label and detail
unitsToPlot = [133 34 119];

figure
rewColors = getDistinctColors('SELECT_ORDER',8);
lw = 2;
indsToPlot_delay = find(time_delay >= -170 & time_delay <=500);
indsToPlot_GC = time_GC >= -200 & time_GC <= 100;
for u = 1:length(unitsToPlot)
    subplot(length(unitsToPlot),3,3*(u-1)+(1:2)); hold on
    unit = unitsToPlot(u);
    maxVal = -inf; minVal = inf;
    for r = 1:nrewards
        x = time_delay(indsToPlot_delay);
        y = delayTraj_byRew_nanmean{r}(unit,indsToPlot_delay);
        err = delayTraj_byRew_nansem{r}(unit,indsToPlot_delay);
        shadedErrorBar(x,y,err,...
            'lineProps',{'r-','color',rewColors{r},'linewidth',lw})
        maxVal = max(y(:)+err(:),maxVal);
        minVal = min(y(:)-err(:),minVal);
    end; clear r
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    xticks([0 100])
    yticks([])
    axis([-inf inf -inf inf])
    
    subplot(length(unitsToPlot),3,3*u); hold on
    for r = 1:nrewards
        x = time_GC(indsToPlot_GC);
        y = GCTraj_byRew_nanmean{r}(unit,indsToPlot_GC);
        err = GCTraj_byRew_nansem{r}(unit,indsToPlot_GC);
        shadedErrorBar(x,y,err,...
            'lineProps',{'r-','color',rewColors{r},'linewidth',lw})
        maxVal = max([y(:)+err(:) ; maxVal]);
        minVal = min([y(:)-err(:) ; minVal]);
    end; clear r
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    xticks([0])
    yticks([])
    
    % Round axes to nearest 5 and set them
    maxVal = ceil(maxVal/5)*5+2;
    minVal = floor(minVal/5)*5-2;
    axis([-inf inf minVal maxVal])
    set(gca,'visible','off')
    plot(0,minVal,'sg','linewidth',1,'markeredgecolor',[0 1 0]*0.8,'markerfacecolor',[0 1 0])

    subplot(length(unitsToPlot),3,3*(u-1)+(1:2)); hold on
    axis([-inf inf minVal maxVal])
    set(gca,'visible','off')
    plot(0,minVal,'sm','linewidth',1,'markeredgecolor',[1 0 1]*0.8,'markerfacecolor',[1 0 1])
    if u == length(unitsToPlot)
        barLength = 100;
        plot(time_delay(indsToPlot_delay(1:barLength/binSize)),minVal*ones(barLength/binSize,1),'k-','linewidth',2)
    end
    plot(min(time_delay(indsToPlot_delay))*ones(100,1),linspace(15,20,100),'k-','linewidth',2)
end; clear u
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 353 443])


% Save it!
figname = ['Fig1C_PSTHExamples_unlabeled']
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);



%% Show a 3D plot of these
load('NeuralDataForFigs_Earl_20220826_152605_allDelays')
rewMethod = 'PCA';
curData = singleNeuronData_pGC(:,unitsToPlot);
curData = curData(:,[2 3 1]); % Show 1 on z, 2 on x, 3 on y
curData_byRew_mean = cellfun(@(x) mean(x), groupDataByLabel(curData,rewardLabels), 'UniformOutput', false);
curData_byRew_mean = reshape([curData_byRew_mean{:}],[3 4])'; % rew x neuron

% Identify the reward axis just in this data
[wR,muR,rewProjData,eigVlsR] = getRewardAxis(curData,rewardLabels,directionLabels,rewMethod);

figure; hold on
set(gcf,'renderer','painters')

scale = 7.5; % how much to stretch the unit vector
ms = 20;
lw = 1.5;
alphas = 0.25*[0.3 0.5 0.3 0.5];
rewColors = getDistinctColors('SELECT_ORDER',8);
% for r = 1:nrewards
%     pairVals = [curData_byRew_mean(r,:); curData_byRew_mean(r,:)];
%     plot3([1000 ; pairVals(2,1)],pairVals(:,2),pairVals(:,3),'k-','color',rewColors{r})
%     plot3(pairVals(:,1),[1000 ; pairVals(2,2)],pairVals(:,3),'k-','color',rewColors{r})
%     plot3(pairVals(:,1),pairVals(:,2),[0 ; pairVals(2,3)],'k-','color',rewColors{r})
% end; clear r
% for r = 1:nrewards
%     curInds = rewardLabels==r;
%     x = curData(curInds,1);
%     y = curData(curInds,2);
%     z = curData(curInds,3);
%     scatter3(x,y,z,ms1,rewColors{r},'filled','markerfacealpha',alphas(r),'markeredgealpha',alphas(r))
% end; clear r

a = quiver3(muR(1),muR(2),muR(3),scale*wR(1),scale*wR(2),scale*wR(3),'color',[0 0.8 0],...
    'linewidth',lw,'maxheadsize',1);
a = quiver3(muR(1),muR(2),muR(3),-scale*wR(1),-scale*wR(2),-scale*wR(3),'color',[0 0.8 0],...
    'linewidth',lw,'maxheadsize',1);
for r = 1:nrewards
    plot3(curData_byRew_mean(r,1),curData_byRew_mean(r,2),curData_byRew_mean(r,3),'o',...
        'color',rewColors{r},'markerfacecolor',rewColors{r},'markersize',ms)
end; clear r
grid on
view(-40,55)
axis([5 20 15 30 15 30])
% xticks([5 20])
% yticks([15 30])
% zticks([15 30])
set(gca,'fontsize',14,'fontname','arial','tickdir','out');
xlabel('Neuron 2')
ylabel('Neuron 3')
zlabel('Neuron 1 (spikes/s)')

% Save it!
figname = ['Fig1D_examplesStateSpace_large']
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Try again, but smaller
figure; hold on
set(gcf,'renderer','painters')

scale = 12; % how much to stretch the unit vector
ms = 12.5;
lw = 2;
alphas = 0.25*[0.3 0.5 0.3 0.5];
rewColors = getDistinctColors('SELECT_ORDER',8);

for r = 1:nrewards
    plot3(curData_byRew_mean(r,1),curData_byRew_mean(r,2),curData_byRew_mean(r,3),'o',...
        'color',rewColors{r},'markerfacecolor',rewColors{r},'markersize',ms)
end; clear r
pts = [muR-scale*wR'; muR+scale*wR'];
plot3(pts(:,1),pts(:,2),pts(:,3),'-','color',[0 0.8 0],'linewidth',lw);
grid on
view(-45,10)
xvals = [2.5 22.5];
yvals = [12.5 32.5];
zvals = [12.5 32.5];
axis([xvals yvals zvals])
xticks(xvals)
yticks(yvals)
zticks(zvals)
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'fontsize',12,'fontname','arial','tickdir','out');
xlabel('Neuron 2')
ylabel('Neuron 3')
zlabel('Neuron 1')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 258 216])

% Save it!
figname = ['Fig1D_exampleStateSpace']
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);