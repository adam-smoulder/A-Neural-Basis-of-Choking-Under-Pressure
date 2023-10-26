%%% This script looks at the angle between the target plane and reward axis
%%% and compares that to a null distribution drawn from the covariance of
%%% the data.
%%%
%%% Adam Smoulder, 7/26/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220713_144001';
    'NeuralDataForFigs_Prez_20220713_144024';
    'NeuralDataForFigs_Rocky_20220713_145614';
    };
nsubjects = length(filename_bySubject);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;
rng('default'); rng(3195)

%% For each animal, get the target plane
rewMethod = 'PCA'
targMethod = 'PCA'
subjectNames = cell(nsubjects,1);
nfactors_bySubj = nan(nsubjects,1);
varExpDir_subjSplit = cell(nsubjects,1);
targProj_subjSplit = cell(nsubjects,1);
rewTargProj_subjSplit = cell(nsubjects,1);
directionLabels_subjSplit = cell(nsubjects,1);
rewardLabels_subjSplit = cell(nsubjects,1);
wRwDAngle_bySubj = nan(nsubjects,1);
nnull = 50000;
nnullalign = 1000;
nullAngles_bySubj = nan(nsubjects,nnull);
nullAlignedAngles_bySubj = nan(nsubjects,nnullalign);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    nfactors = size(binnedData_pGC,2);
    nfactors_bySubj(f) = nfactors;

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
    
    % Get target plane and projection
    [wD,muD,~,varExpDir_subjSplit{f}] = getTargetPlane(binnedData_pGC(~skipLabels,:),...
        rewardLabels(~skipLabels,:),directionLabels(~skipLabels,:),targMethod);
    targProj_subjSplit{f} = (binnedData_pGC-muD)*wD;
    directionLabels_subjSplit{f} = directionLabels;
    rewardLabels_subjSplit{f} = rewardLabels;
    
    % Show angle between top 2 dims of wD
    disp(['Angle between target plane dims (deg): ' num2str(acosd(wD(:,1)'*wD(:,2)))])
    
    % Get reward axis
    [wR,muR] = getRewardAxis(binnedData_pGC(~skipLabels,:),...
        rewardLabels(~skipLabels,:),directionLabels(~skipLabels,:),rewMethod);
    
    % Orthogonalize the reward axis w.r.t. the target plane
    [wDR,~] = qr([wD wR]);
    wDR = wDR(:,1:3);
    wDR(:,1:2) = wD; % assume wD was already orthogonal
    
    % Get the angle between the original reward axis and the target plane
    rewAlterAngle = 180/pi*allProjectionAngles(wDR'*wR,[0 0 1]'); % Angle between rew axis and the orthogonal rew axis
    wRwDAngle_bySubj(f) = 90-min(rewAlterAngle,180-rewAlterAngle);
   
    % Get the null distribution for these using the Elsayed paper method
    sigma = cov(binnedData_pGC);
    [U,S,~] = svd(sigma);
    nullAngles = nan(nnull,1);
    for i = 1:nnull
        curVec = randn([nfactors 1]);
        [curVec,~,~] = svd(U*sqrt(S)*curVec/norm(U*sqrt(S)*curVec));
        nullAngles(i) = subspace(wD,curVec(:,1))*180/pi; % assumes both are normalized already!!!
    end; clear i
    nullAngles_bySubj(f,:) = nullAngles;
    
    % To get a null distribution of aligned vectors, we'll randomly split
    % the data in two, calculate a target plane in 1, and calculate a 1st
    % target plane dim in the other
    dirLabels = directionLabels(~skipLabels);
    rewLabels = rewardLabels(~skipLabels);
    data = binnedData_pGC(~skipLabels,:);
    [~,inds_byDir_byRew,n_byDir_byRew] = groupDataByLabel(dirLabels,[dirLabels rewLabels]);
    for b = 1:nnullalign
        repInds = cellfun(@(x) x(randperm(length(x))), inds_byDir_byRew, 'UniformOutput', false);
        v1Inds = cellArrayToVector(cellfun(@(x) x(1:floor(length(x)/2)), repInds,'uniformoutput',false));
        v2Inds = cellArrayToVector(cellfun(@(x) x(floor(length(x)/2)+1:end), repInds,'uniformoutput',false));
        wDv1 = getTargetPlane(data(v1Inds,:),rewLabels(v1Inds,:),dirLabels(v1Inds,:),rewMethod);
        wDv2 = getTargetPlane(data(v2Inds,:),rewLabels(v2Inds,:),dirLabels(v2Inds,:),rewMethod);
        nullAlignedAngles_bySubj(f,b) = abs(subspace(wDv1,wDv2(:,1))*180/pi);
    end; clear b
    
%     dirLabels = directionLabels(~skipLabels);
%     rewLabels = rewardLabels(~skipLabels);
%     data = binnedData_pGC(~skipLabels,:);
%     [~,inds_byDir_byRew,n_byDir_byRew] = groupDataByLabel(dirLabels,[dirLabels rewLabels]);
%     for b = 1:nnullalign
%         repInds = cellfun(@(x) x(randperm(length(x))), inds_byDir_byRew, 'UniformOutput', false);
%         v1Inds = cellArrayToVector(cellfun(@(x) x(1:floor(length(x)/2)), repInds,'uniformoutput',false));
%         v2Inds = cellArrayToVector(cellfun(@(x) x(floor(length(x)/2)+1:end), repInds,'uniformoutput',false));
%         wRv1 = getRewardAxis(data(v1Inds,:),rewLabels(v1Inds,:),dirLabels(v1Inds,:),rewMethod);
%         wRv2 = getRewardAxis(data(v2Inds,:),rewLabels(v2Inds,:),dirLabels(v2Inds,:),rewMethod);
%         nullAlignedAngles_bySubj(f,b) = abs(subspace(wRv1,wRv2)*180/pi);
%     end; clear b
end; clear f
disp('Done processing data')



%% Show the null distributions and reward axis values for each animal
figure
binEdges = linspace(0,90,40); % 0.5deg bins
for f = 1:nsubjects
    % Plot the null distribution
    subplot(1,nsubjects,f); hold on
    h = histogram(nullAngles_bySubj(f,:),binEdges,'normalization','probability');
    h.FaceColor = [0.3 0.3 0.3];
    h2 = histogram(nullAlignedAngles_bySubj(f,:),binEdges,'normalization','probability');
    h2.FaceColor = [0.9 0.9 0.9];
    set(gca,'fontsize',14,'fontname','arial','tickdir','out');
    title(['Monkey ' subjectNames{f}(1) ' (' num2str(nfactors_bySubj(f)) ' factors)'])
    xticks([0 45 90])
%     yticks([0 0.05 0.10 0.15])
    axis([0 90 0 0.85])
    if f == 1
        ylabel('Relative Frequency')
    else
        yticklabels([])
    end
    xlabel('Angle (degrees)')
    
    % Plot the actual value
    plot(wRwDAngle_bySubj(f)*ones(100,1),linspace(0,0.78,100),'k-','linewidth',1)
    curPrctile = round(sum(wRwDAngle_bySubj(f)>nullAngles_bySubj(f,:))/nnull*100);
    text(wRwDAngle_bySubj(f),0.84,[num2str(round(wRwDAngle_bySubj(f))) '^o'],...
        'horizontalalignment','center','fontname','arial','fontsize',12)
    text(wRwDAngle_bySubj(f),0.8,['(' num2str(curPrctile) ' prctile of null dist)'],...
        'horizontalalignment','center','fontname','arial','fontsize',12)
    
    if f == 1
        legend('Random','Same','Data','location','nw')
    end
end; clear f

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 1642 420])

% Save it
figname = ['FigS2_angleBetweenTargetPlaneAndRewAxis_rew' rewMethod '_targ' targMethod ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Make chibi plot for inclusion in Fig 2

figure
step = 5;
binEdges = -step/2:step:(90+step); % 0.5deg bins
x = 0.5*(binEdges(1:end-1) + binEdges(2:end));
for f = 1:nsubjects
    % Plot the null distribution
    subplot(nsubjects,1,f); hold on
%     h = histogram(nullAngles_bySubj(f,:),binEdges,'normalization','probability');
%     h.FaceColor = [0.3 0.3 0.3];
    h = histcounts(nullAngles_bySubj(f,:),binEdges,'normalization','probability');
    plot(x,h,'-','linewidth',1.5,'color',[0.2 0.2 0.2])
%     h2 = histogram(nullAlignedAngles_bySubj(f,:),binEdges,'normalization','probability');
%     h2.FaceColor = [0.9 0.9 0.9];
    set(gca,'fontsize',9,'fontname','arial','tickdir','out');
    xticks([0 90])
    axis([0 90 0 0.35])
    ylabel('Frequency')
    yticklabels([])
    xlabel('Angle')
    
    % Plot the actual value
    xline(wRwDAngle_bySubj(f),'k-','linewidth',2.5) % for some odd reason it's transparent-ed?
    xline(wRwDAngle_bySubj(f),'k-','linewidth',2.5)
    xline(wRwDAngle_bySubj(f),'k-','linewidth',2.5)
    
%     if f == 1
%         legend({'Null','Data'},'location','nw')
%     end
    
    curPrctile = round(sum(wRwDAngle_bySubj(f)>nullAngles_bySubj(f,:))/nnull*100);
    disp([num2str(round(wRwDAngle_bySubj(f)))])
    disp(['(' num2str(curPrctile) ' prctile of null dist)'])
end; clear f

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 131 372])

% Save it
figname = ['FigS2_MINIangleBetweenTargetPlaneAndRewAxis_rew' rewMethod '_targ' targMethod ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);


%% Try another version that's between the two
vals_bySubj = [0.6, 0.4, 0.9];

figure
binEdges = linspace(0,90,40); % 0.5deg bins
for f = 1:nsubjects
    val = vals_bySubj(f);
    
    % Plot the null distribution
    subplot(nsubjects,1,f); hold on
    h = histogram(nullAngles_bySubj(f,:),binEdges,'normalization','probability');
    h.FaceColor = [0.3 0.3 0.3];
    h2 = histogram(nullAlignedAngles_bySubj(f,:),binEdges,'normalization','probability');
    h2.FaceColor = [0.9 0.9 0.9];
    set(gca,'fontsize',12,'fontname','arial','tickdir','out');
    title(['Monkey ' subjectNames{f}(1) ' (' num2str(nfactors_bySubj(f)) ' factors)'])
    xticks([0 90])
%     yticks([0 0.05 0.10 0.15])
    axis([0 90 0 val])
%     if f == 1
    ylabel('Relative Frequency')
%     else
%         yticklabels([])
%     end
    xlabel('Angle')
    
    % Plot the actual value
    plot(wRwDAngle_bySubj(f)*ones(100,1),linspace(0,val*0.8,100),'k-','linewidth',1)
    curPrctile = round(sum(wRwDAngle_bySubj(f)>nullAngles_bySubj(f,:))/nnull*100);
    text(wRwDAngle_bySubj(f),val*0.9647,[num2str(round(wRwDAngle_bySubj(f))) '^o'],...
        'horizontalalignment','center','fontname','arial','fontsize',12)
    text(wRwDAngle_bySubj(f),val*0.8824,['(' num2str(curPrctile) '%ile of null)'],...
        'horizontalalignment','center','fontname','arial','fontsize',12)
    
%     if f == 1
%         legend('Random','Same','Data','location','nw')
%     end
end; clear f

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[0 500] 359 950])

% Save it
figname = ['FigS2_MIDangleBetweenTargetPlaneAndRewAxis_rew' rewMethod '_targ' targMethod ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);
