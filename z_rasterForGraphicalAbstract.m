addpath(genpath('C:\Users\Smoulgari\Documents\MATLAB\'))
load('Rocky_day01_half1_stitchedResults_14-Apr-2022085224')


%%
window = -600:400;
i = 111;
spikeMat = dayTrialData(i).neuralData.spikeMatrix;
time = dayTrialData(i).time;
GCInd = find(time==dayTrialData(i).t3_goCueTime);

curSpikeMat = spikeMat(1:150,GCInd + window); % Stick to anterior array; good delta

% Pick the 25 units with most spikes and randomize their order
nunitsToKeep = 40;
[counts,order] = sort(sum(curSpikeMat(:,end-200:end),2),'descend');
order = order(1:nunitsToKeep);
order = order(randperm(nunitsToKeep));
curSpikeMat = curSpikeMat(order,:);

figure; hold on
set(gcf,'position',[745 458 426 268])
xline(find(window==-150),'r-')
xline(find(window==50),'r-')
xline(find(window==0),'g-')
set(gca,'visible','off')

imagesc(1-curSpikeMat)
axis([1 length(window) 1 nunitsToKeep])
colormap gray

figfolder = 'C:\Users\Smoulgari\Documents\MATLAB\~chokingUnderPressure\~forManuscript\';
figname = 'rasterForGraphicalAbstract'
saveFigAndSvg([figfolder 'recentFigs/'],figname)
saveFigAndSvg([figfolder 'allFigs/'],[figname '_' grabDateTimeString])
