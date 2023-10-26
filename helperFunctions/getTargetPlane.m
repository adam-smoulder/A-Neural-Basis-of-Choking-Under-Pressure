function [wD,muD,targProjData,modelQuant] = getTargetPlane(data,rewardLabels,directionLabels,method,varargin)
% This function identifies the target plane from input data and direction
% labels after balancing by reward.
%
% Inputs:
% - data: [ntrials x nfactors] single neural data bin per trial
% - rewardLabels: [ntrials x 1] categorical reward cue label per trial
% - directionLabels: [ntrials x 1] categorical direction label per trial
% - method: string, decides which method to calculate the target plane with
% Options include:
%   - LR: Use linear regression on dir/rew balanced data against target
%   locations
%       DETAILS HERE
%   - LDA: Use LDA on dir/rew-balanced data against direction labels
%       DETAILS HERE
%   - PCA: Use PCA on reward-balanced direction means
%       DETAILS HERE
%
% Optional inputs:
% - ndimout: scalar, number of dimensions for output to have. Only applies
% to PCA/LDA methods. Default value is 2.
% - orthogonalize: boolean, do we orthogonalize wD using QR decomposition? 
% Default value is true.
%
% Outputs:
% - wD: [nfactors x numdimout] target plane
% - muD: [1 x nfactors] mean of data used for model fitting
% - targProjData: [ntrials x numdimout] projection of single trial data
% onto the target plane
%   - For LR, we include the offset in this
% - modelQuant: value(s) reflecting a metric that differs depending on the 
% model used for fitting:
%   - PCA methods: percent variance explained by each dimension
%   - LDA methods: the model vector's corresponding eigenvalues (SNR)
%   - LR methods: R^2 (not cross-validated; only to be used for reference)
%
% Adam Smoulder, 7/19/22

% Assign optional arguments
if isempty(varargin)
    ndimout = 2;
    orthogonalize = true;
elseif length(varargin)==1
    ndimout = varargin{1};
    orthogonalize = true;
else
    ndimout = varargin{1};
    orthogonalize = varargin{2};
end

% Group data by reward and direction
[ntrials,nfactors] = size(data);
[data_byDir_byRew,inds_byDir_byRew] = groupDataByLabel(data,[directionLabels,rewardLabels]);
[ndirections,nrewards] = size(data_byDir_byRew);

% For LDA/LR, we'll need to oversample to balance dir/rew conditions
if contains(method,'LDA') || contains(method,'LR')
    if size(data,2) < 100
        n_OS = 5000; % # points to oversample per cond. Bigger # = takes longer but more stable
    else
        n_OS = 1000;
    end
%     OSInds = cellArrayToVector(cellfun(@(x) x(randi(length(x),[n_OS 1])),...
%         inds_byDir_byRew,'uniformoutput',false)); % true random
    OSInds = cellArrayToVector(cellfun(@(x) x(permutationOversample(length(x),n_OS)),...
        inds_byDir_byRew,'uniformoutput',false)); % repeated random
    muD = mean(data(OSInds,:));
end

% From here, fitting depends on method
if contains(method,'PCA')
    % Get means for directions as the avg. over dirs of dir x rew means
    data_byDir_byRew_mean = cellfun(@(x) mean(x,1), data_byDir_byRew,'uniformoutput',false);
    for d = 1:ndirections, for r = 1:nrewards
            if isempty(data_byDir_byRew_mean{d,r}), data_byDir_byRew_mean{d,r} = nan(1,nfactors); end
    end; clear r; end; clear d
    data_byDir_byRew_mean = ... % go from cell to matrix
        reshape([data_byDir_byRew_mean{:}],[nfactors ndirections size(data_byDir_byRew,2)]);
    data_byDir_mean = nanmean(data_byDir_byRew_mean,3)'; % nrewards x nfactors
    
    % Do PCA to get the target plane
    muD = mean(data_byDir_mean);
    [wD,z,~,~,modelQuant] = pca(data_byDir_mean-muD,'numcomponents',ndimout); % z is projection of direction means for debugging
    LROffset = zeros(1,ndimout);
    
elseif contains(method,'LDA')
    % Fit LDA model with oversampled input
    [wD,~,modelQuant,cvaccs] = linearDiscAnalysis(data(OSInds,:),directionLabels(OSInds),'NUMDIMOUT',ndimout); % cvaccs is cv'd decode acc for debugging
    LROffset = zeros(1,ndimout);
    
elseif contains(method,'LR')
    % Set target locations with radius 1
    theta = 0:45:315;
    targetLocs = [cosd(theta)' sind(theta)'];
    
    % Fit LR model with oversampled input
    regInput = [data(OSInds,:)-muD ones(length(OSInds),1)]; % add ones for offset
    [wD,~,resid] = mvregress(regInput,targetLocs(directionLabels(OSInds),:));
    LROffset = wD(end,:);
    wD = wD(1:end-1,:); % remove offset term; we don't need it
    modelQuant = 1-var(resid)./var(targetLocs(directionLabels(OSInds),:)); % R^2
    
else
    error('What method are you using?')
end

% Orthogonalize wD if desired
wD = wD(:,1:ndimout);
if orthogonalize
    [wD,~] = qr(wD);
    wD = wD(:,1:ndimout);
end

% Project single trial data onto the target plane
targProjData = (data-muD)*wD;

% Orient wD such that up-right targets are pointed towards (1,1), at least
% for the top 2 dims
if ~(contains(method,'LR') && ~orthogonalize) % unorthogonalized LR models are already oriented
    % First reflect if needed (if up-right is CW of down-right)
    upRightProjMean = mean(targProjData(directionLabels==2,1:2));
    upRightProjUnitVec = upRightProjMean/norm(upRightProjMean);
    downRightProjMean = mean(targProjData(directionLabels==8,1:2));
    downRightProjUnitVec = downRightProjMean/norm(downRightProjMean);
    if downRightProjUnitVec(1)*upRightProjUnitVec(2) < downRightProjUnitVec(2)*upRightProjUnitVec(1)
        upRightProjUnitVec = upRightProjUnitVec.*[1 -1];
        wD(:,1:2) = wD(:,1:2).*[1 -1];
        targProjData(:,1:2) = targProjData(:,1:2).*[1 -1];
    end
    
    % Now get the rotation angle
    newUnitVec = 1./sqrt([2 2]);
    rotAngle = acosd(upRightProjUnitVec*newUnitVec');
    rotMat = [cosd(rotAngle) sind(rotAngle) ; -sind(rotAngle) cosd(rotAngle)]; % cur data * this = pointed to new vec
    if ~(norm(upRightProjUnitVec*rotMat-newUnitVec) < 1E-8) % wrong angle directon; transpose
        rotMat = rotMat';
    end
    assert((norm(upRightProjUnitVec*rotMat-newUnitVec) < 1E-8));
    
    % Rotate the projection
    wD(:,1:2) = wD(:,1:2)*rotMat;
    targProjData(:,1:2) = targProjData(:,1:2)*rotMat;
else % If LR and no orthogonalizing, add the offset
    targProjData = targProjData+LROffset;
end


end

