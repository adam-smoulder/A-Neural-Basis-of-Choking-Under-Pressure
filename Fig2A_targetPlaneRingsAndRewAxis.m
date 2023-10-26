%%% This script makes target plane rings with axes defined by PCA.
%%%
%%% Adam Smoulder, 7/19/22

filename_bySubject = {...
    'NeuralDataForFigs_Earl_20220826_152605_allDelays';
    'NeuralDataForFigs_Prez_20220826_152628_allDelays';
    'NeuralDataForFigs_Rocky_20220826_152918_allDelays';
    };
nsubjects = length(filename_bySubject);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;

%% For each animal, get the target plane
rewMethod = 'PCA'
targMethod = 'PCA'
subjectNames = cell(nsubjects,1);
dirPve_subjSplit = cell(nsubjects,1);
targProj_subjSplit = cell(nsubjects,1);
rewTargProj_subjSplit = cell(nsubjects,1);
directionLabels_subjSplit = cell(nsubjects,1);
rewardLabels_subjSplit = cell(nsubjects,1);
wRwDAngle_bySubj = nan(nsubjects,1);
projPVE_subjSplit = cell(nsubjects,1);
for f = 1:nsubjects
    % Load data
    disp(['Analyzing data for subject ' num2str(f)])
    load(filename_bySubject{f})
    subjectNames{f} = taskInfo.subjectName;
    nfactors = size(binnedData_pGC,2);

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
    [wD,muD,~,dirPve_subjSplit{f}] = getTargetPlane(binnedData_pGC(~skipLabels,:),...
        rewardLabels(~skipLabels),directionLabels(~skipLabels),targMethod);
    targProj_subjSplit{f} = (binnedData_pGC(~skipLabels,:)-muD)*wD;
    directionLabels_subjSplit{f} = directionLabels(~skipLabels);
    rewardLabels_subjSplit{f} = rewardLabels(~skipLabels);
    
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

    % Assert the direction of the reward axis again
    tempRewProj = (binnedData_pGC-muR)*wDR(:,3);
    if mean(tempRewProj(rewardLabels==1)) > mean(tempRewProj(rewardLabels==3))
        wDR(:,3) = -wDR(:,3);
    end
    
    % Store it
    rewTargProj_subjSplit{f} = [(binnedData_pGC(~skipLabels,:)-muD)*wDR(:,1:2) (binnedData_pGC(~skipLabels,:)-muR)*wDR(:,3)];
    
    % Calculate the variance explained by each dimension
    data_byDir_byRew = groupDataByLabel(binnedData_pGC(~skipLabels,:),[directionLabels(~skipLabels) rewardLabels(~skipLabels)]);
    dataMeans_byDir_byRew = cellfun(@(x) mean(x), data_byDir_byRew, 'UniformOutput', false);
    dataMeans_byDir_byRew = reshape(cellArrayToVector(dataMeans_byDir_byRew),[nfactors ndirections nrewards]);
    dataMeans_byCond = dataMeans_byDir_byRew(:,:)'; % (ndir*nrew) x nfactors
    totalVar = sum(var(dataMeans_byCond));
    projDataMeans_byCond = [(dataMeans_byCond-muD)*wDR(:,1:2) (dataMeans_byCond-muR)*wDR(:,3)];
    projPVE_subjSplit{f} = var(projDataMeans_byCond)./totalVar*100;
end; clear f
disp('Done processing data')


%% We also want some smaller plots for insets to emphasize the rings inverted-U

figure
set(gcf,'renderer','painters')
rewColors = getDistinctColors('SELECT_ORDER',8);
dirColors = getDistinctColors('SELECT_ORDER',16);
bgColor = [0.9686 0.9059 0.8078]+0.03; % used for inset
% bgColor = [0 0 0];
ms1 = 9;
ms2 = 4;
lw1 = 3;
lw2 = 3;
insetAxisValues_subjSplit = {...
    [-58 -22 -58 -22]; % E-green
    [30 65 25 60]; % P-plum
    [-143 -107 13 49]; % R-green
    };
subplotOrder_subjSplit = [...
    {{[51 52 [51 52]+12],[1:4 (1:4)+12 (1:4)+24 (1:4)+36]}},...
    {{[55 56 [55 56]+12],[(1:4)+4 (1:4)+16 (1:4)+28 (1:4)+40]}},...
    {{[59 60 [59 60]+12],[(1:4)+8 (1:4)+20 (1:4)+32 (1:4)+44]}}
    ];



% Transform to the new subplot orders
W = 12;
L = 6;
subplotOrder_subjSplit = cellfun(@(y) cellfun(@(x) ...
    sort(flat((2*x-1+floor((x-0.1)/W)*2*W)+[0 1 2*W 2*W+1]'))',y,'uniformoutput',false),...
    subplotOrder_subjSplit,'uniformoutput',false);
for f = 1:nsubjects
    subplotOrder_subjSplit{f}{1}(end-3:end) = [];
end; clear f
if strcmp(targMethod,'PCA')
    view_subjSplit = {[65 10],[-140 10],[-190 17]}; % good for PCA axis equal (orig)
    tickVals_subjSplit = {[-210:35:210],[-210:35:210],[-240:60:240]};
    axis_subjSplit = {[-70 105 -105 70 -70 70],[-75 75 -100 100 -70 105],[-180 180 -240 180 -75 65]};
else
    view_subjSplit = {[50 10],[-150 28],[-70 25]}; % good for PCA axis equal
    tickVals_subjSplit = {[-210:35:210],[-210:35:210],[-240:60:240]};
    axis_subjSplit = {[-70 70 -70 70 -70 70],[-75 75 -105 105 -70 105],[-180 180 -200 200 -140 140]};
end

outerPos_subjSplit = {[0.09,0.4,0.29,0.5],[0.40,0.38,0.26,0.59],[0.635,0.42,0.32,0.46]};



for f = 1:nsubjects
    % Get the current points
    targProj_byDir_byRew = groupDataByLabel(targProj_subjSplit{f},[directionLabels_subjSplit{f} rewardLabels_subjSplit{f}]);
    targProj_byDir_byRew_means = cellfun(@(x) mean(x), targProj_byDir_byRew,'uniformoutput',false);
    [ndirections,nrewards] = size(targProj_byDir_byRew_means);
    directions = unique(directionLabels_subjSplit{f}); % diff for Prez
    temp = reshape([targProj_byDir_byRew_means{:}],[2,ndirections,nrewards]);
    axisMax = 1.2*max(temp(:,:),[],2);
    axisMin = 1.2*min(temp(:,:),[],2);
    
    % Plot insets
    fullMax = max(abs([axisMin ; axisMax]))/1.1;
    subplot(12,24,subplotOrder_subjSplit{f}{1}); hold on
    z = axis_subjSplit{f}(5)+1;
    P = [insetAxisValues_subjSplit{f}([1 3]) z; ...
        insetAxisValues_subjSplit{f}([1 4]) z;...
        insetAxisValues_subjSplit{f}([2 4]) z;...
        insetAxisValues_subjSplit{f}([2 3]) z;...
        ];
    patch(P(:,1),P(:,2),P(:,3),bgColor);
    for r = 1:nrewards
        curMeans = reshape([targProj_byDir_byRew_means{:,r}],[2 ndirections])';
        curMeans = curMeans-mean(curMeans); % we can center within rewards when projecting in target plane
        ringPoints = [curMeans ; curMeans(1,:)];
        plot(ringPoints(:,1),ringPoints(:,2),'-','linewidth',1.5,...
            'color',rewColors{r})
        for d = 1:ndirections
            plot(curMeans(d,1),curMeans(d,2),'o','markersize',1.5*ms1,'linewidth',3.5,...
                'markerfacecolor',dirColors{directions(d)},'markeredgecolor',rewColors{rewards(r)});
        end; clear d
    end; clear r
    axis(insetAxisValues_subjSplit{f});
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    xlabel('Target Axis 1')
    if f == 1
        ylabel('Target Axis 2')
    end
    xticks([])
    yticks([])
    
    
    
    
    % Get current points
    subplot(12,24,subplotOrder_subjSplit{f}{end}); hold on
    rewTargProj_byDir_byRew = groupDataByLabel(rewTargProj_subjSplit{f},[directionLabels_subjSplit{f} rewardLabels_subjSplit{f}]);
    rewTargProj_byDir_byRew_means = cellfun(@(x) mean(x), rewTargProj_byDir_byRew,'uniformoutput',false);
    [ndirections,nrewards] = size(rewTargProj_byDir_byRew_means);
    directions = unique(directionLabels_subjSplit{f}); % diff for Prez
    for r = 1:nrewards
        curMeans = reshape([rewTargProj_byDir_byRew_means{:,r}],[3 ndirections])';
        curMeans(:,1:2) = curMeans(:,1:2)-mean(curMeans(:,1:2)); % we can center within rewards when projecting in target plane
        ringPoints = [curMeans ; curMeans(1,:)];
        plot3(ringPoints(:,1),ringPoints(:,2),ringPoints(:,3),'-','linewidth',1.5*lw1,...
            'color',rewColors{r})
        for d = 1:ndirections
%             plot3(curMeans(d,1),curMeans(d,2),curMeans(d,3),'o','markersize',1.25*ms1,'linewidth',lw1,...
%                 'markeredgecolor',dirColors{directions(d)},'markerfacecolor',dirColors{directions(d)});
            plot3(curMeans(d,1),curMeans(d,2),curMeans(d,3),'o','markersize',1.25*ms1,'linewidth',lw2,...
                'markeredgecolor',rewColors{rewards(r)},'markerfacecolor',dirColors{directions(d)});
        end; clear d
    end; clear r
    set(gca,'fontsize',14,'fontname','arial','tickdir','out');
    if f == 1
        zlabel('Reward Axis') 
    end
    xlabel('Targ. Ax. 1')
    ylabel('Targ. Ax. 2')
    title(['Monkey ' subjectNames{f}(1)])
    view(view_subjSplit{f})
    if f==2
        axis equal
    end
    axis(axis_subjSplit{f})
    set(gca,'outerposition',outerPos_subjSplit{f})

% %     r = rectangle('Position',rectangleVals,'facecolor',bgColor,'edgecolor','k','linewidth',2)
% 
%     % Plot the prism indicating the inset
% %     z1 = axis_subjSplit{f}(5)+1;
% %     z2 = axis_subjSplit{f}(6)+10;
% %     P = [insetAxisValues_subjSplit{f}([1 3]) z1; ...
% %         insetAxisValues_subjSplit{f}([1 4]) z1;...
% %         insetAxisValues_subjSplit{f}([2 4]) z1;...
% %         insetAxisValues_subjSplit{f}([2 3]) z1;...
% %         insetAxisValues_subjSplit{f}([1 3]) z1;...
% %         insetAxisValues_subjSplit{f}([1 3]) z2; ...
% %         insetAxisValues_subjSplit{f}([1 4]) z2; ...
% %         insetAxisValues_subjSplit{f}([1 4]) z1; ...
% %         insetAxisValues_subjSplit{f}([1 4]) z2; ...
% %         insetAxisValues_subjSplit{f}([2 4]) z2;...
% %         insetAxisValues_subjSplit{f}([2 4]) z1;...
% %         insetAxisValues_subjSplit{f}([2 4]) z2;...
% %         insetAxisValues_subjSplit{f}([2 3]) z2;...
% %         insetAxisValues_subjSplit{f}([2 3]) z1;...
% %         insetAxisValues_subjSplit{f}([2 3]) z2;...
% %         insetAxisValues_subjSplit{f}([1 3]) z2;...
% %         ];
% %     p = patch(P(:,1),P(:,2),P(:,3),bgColor,'edgecolor','k','linewidth',0.5,'faceAlpha',1);
%     p = patch(P(:,1),P(:,2),P(:,3),bgColor,'edgecolor','k','linewidth',0.5,'faceAlpha',1);
%     
%     
%     z1 = axis_subjSplit{f}(5)+1;
%     z2 = axis_subjSplit{f}(6)+10;
%     
%     plot3(mean(insetAxisValues_subjSplit{f}([1 2]))*ones(100,1),...
%         mean(insetAxisValues_subjSplit{f}([3 4]))*ones(100,1),...
%         linspace(z1,z2,100),'-','color','k','linewidth',1.5)
    grid on
    xticklabels([])
    yticklabels([])
    zticklabels([])
    xticks(tickVals_subjSplit{f})
    yticks(tickVals_subjSplit{f})
    zticks(tickVals_subjSplit{f})
    
    if f == 3 % Draw dotted line
        
    end; clear f
    
end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)-[500 500] 1680 933])

% Save it
figname = ['Fig2A_targetPlaneRingsAndRewAxis_withInsets_targ' targMethod '_rew' rewMethod ]
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);