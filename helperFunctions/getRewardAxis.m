function [wR,muR,rewProjData,modelQuant] = getRewardAxis(data,rewardLabels,directionLabels,method,varargin)
% This function identifies the reward axis from input data and reward
% labels after balancing by direction.
%
% Inputs:
% - data: [ntrials x nfactors] single neural data bin per trial
% - rewardLabels: [ntrials x 1] categorical reward cue label per trial
% - directionLabels: [ntrials x 1] categorical direction label per trial
% - method: string, decides which method to use to calculate the reward
% axis. Options include:
%   - 'PCA': Use PCA on the direction-balanced reward means
%       This dimension will be the one which simply captures the most
%       variance of the reward means. Does not try to stifle other sources
%       of variability, but only takes in 4 points.
%   - 'LDA': Use LDA on dir/rew-balanced data against reward labels
%       This dimension will attempt to maximize the distance between the
%       classes by both separating the means (maximizing signal variance)
%       and squashing other sources of variaibility (minimizing noise)
%   - 'LR': Use linear regression on dir/rew balanced data against reward
%   labels
%       This dimension attempts to minimize variability about the
%       categorical reward values of 1-4. Not recommended as it assumes an
%       arbitrary distance of 1 between each reward (which may be
%       nonsensical for Jackpots) and that it assumes order of data
%   - 'PCA_noJ': same as PCA, but fits the model without using Jackpots
%   - 'LDA_noJ': same as LDA, but "     "     "       "       "      "
%   - 'LR_noJ': same as LR, but   "     "     "       "       "      "   
% 
% Optional Inputs:
% - ndimout: scalar, number of dimensions for output to have. Only
% applies to PCA/LDA methods. Default value is 1
% - orthogonalize: boolean, do we orthogonalize wR using QR decomposition?
% Only matters if numdimout > 1 (if 1, we'll always normalize it). Default
% value is true
%
% Outputs:
% - wR: [nfactors x numdimout] reward axis
% - muR: [1 x nfactors] mean of data used for model fitting
% - rewProjData: [ntrials x numdimout] projection of single trial data
% onto the reward axis
% - modelQuant: value(s) reflecting a metric that differs depending on the 
% model used for fitting:
%   - PCA methods: percent variance explained by each dimension
%   - LDA methods: the model vector's corresponding eigenvalues (SNR)
%   - LR methods: R^2 (not cross-validated; only to be used for reference)
%
% Adam Smoulder, 7/19/22 (edit 1/5/23)

% Assign optional arguments
if isempty(varargin)
    ndimout = 1;
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

% If method calls for no jackpots, exclude them here
if contains(method,'noJ')
    disp('Skipping Jackpots for model fitting; assuming J rewards are highest reward label value')
    data_byDir_byRew(:,end) = [];
    inds_byDir_byRew(:,end) = [];
end

% If method only calls for Small and Jackpot, only keep them
if contains(method,'SJ')
    disp('Skipping Medium/Large for model fitting; assuming M/L reward labels 2 and 3')
    data_byDir_byRew(:,[2 3]) = [];
    inds_byDir_byRew(:,[2 3]) = [];
end

% For LDA/LR, we'll need to oversample to balance dir/rew conditions
if contains(method,'LDA') || contains(method,'LR')
    n_OS = 5000; % # points to oversample per cond. Bigger # = takes longer but more stable
    OSInds = cellArrayToVector(cellfun(@(x) x(permutationOversample(length(x),n_OS)),...
        inds_byDir_byRew,'uniformoutput',false)); % repeated random
    muR = mean(data(OSInds,:));
end

% From here, fitting depends on method
if contains(method,'PCA')
    % Get means for rewards as the avg. over dirs of dir x rew means
    data_byDir_byRew_mean = cellfun(@(x) mean(x,1), data_byDir_byRew,'uniformoutput',false);
    data_byDir_byRew_mean = ... % go from cell to matrix
        permute(reshape([data_byDir_byRew_mean{:}],[nfactors ndirections size(data_byDir_byRew,2)]),[1 3 2]);
    data_byRew_mean = mean(data_byDir_byRew_mean,3)'; % nrewards x nfactors
    
    % Do PCA to get the reward axis
    muR = mean(data_byRew_mean);
    [wR,z,~,~,modelQuant] = pca(data_byRew_mean-muR,'numcomponents',ndimout); % z is projection of reward means for debugging
%     [wR,z,modelQuant,~,~] = pca(data_byRew_mean-muR,'numcomponents',ndimout); % z is projection of reward means for debugging

elseif contains(method,'LDA')
    % Fit LDA model with oversampled input
    [wR,~,modelQuant,cvaccs] = linearDiscAnalysis(data(OSInds,:),rewardLabels(OSInds),'NUMDIMOUT',ndimout); % cvaccs is cv'd decode acc for debugging
    
elseif contains(method,'LR')
    % Fit LR model with oversampled input
    regInput = [data(OSInds,:)-muR ones(length(OSInds),1)]; % add ones for offset
    [wR,~,resid] = mvregress(regInput,rewardLabels(OSInds));
    wR = wR(1:end-1,:); % remove offset term; we don't need it
    modelQuant = 1-var(resid)/var(rewardLabels(OSInds)); % R^2
else
    error('What method are you using?')
end

% Orthogonalize (or normalize) wR if desired
wR = wR(:,1:ndimout);
if orthogonalize
    [wR,~] = qr(wR);
    wR = wR(:,1:ndimout);
elseif ndimout == 1
    wR = wR./norm(wR);
end

% Project single trial data onto the reward axis
rewProjData = (data-muR)*wR;

% Orient wR such that Large > Small reward (at least along 1st dim)
if mean(rewProjData(rewardLabels==1),1) > mean(rewProjData(rewardLabels==3),1)
    wR = -wR;
    rewProjData = -rewProjData;
end

end

