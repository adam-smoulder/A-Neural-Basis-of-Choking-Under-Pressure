function [decodeAcc_byRew_mean,decodeAcc_byRew_sem] = nbDecodeHell(binnedData,directionLabels,rewardLabels,nfolds,nreps_mean,nreps_sem)
% This function is a copypasta of what used to be a script that does the
% following:
% - for each reward size individually...
%   - get mean decode accuracy by running many subsamples to whatever the
%   minimum number of trials is through CV Naive Bayes decoding
%   - get SE decode accuracy by running through nested-CVing
%
% Inputs:
% - binnedData: [ntrials x nfactors] data that we want to decode direction
% from.
% - directionLabels: [ntrials x 1] direction label for each trial
% - rewardLabels: [ntrials x 1] reward label for each trial
% - nfolds: scalar, how many folds to use (applies to both)
% - nreps_mean: scalar, how many repeats to use for calculating mean decode
% accuracy. This tends to run quickly, so ~1000 seems typically fair
% - nreps_sem: scalar, how many repeats to use for calculating SE of decode
% accuracy (using nested CV). This is pretty slow, so ~100 is probably
% fair. 
% 
% Outputs:
% - decodeAcc_byRew_mean: [nrewards x 1] average decode accuracy (%)
% - decodeAcc_byRew_sem: [nrewards x 1] SE of decode accuracy (%)
%
% Adam Smoulder, 6/15/21


[~,trialInds_byDir_byRew] = groupDataByLabel(binnedData,[directionLabels,rewardLabels]);
rewards = unique(rewardLabels); nrewards = length(rewardLabels);
directions = unique(directionLabels);
[ndirections,nrewards] = size(trialInds_byDir_byRew);
n_byDir_byRew = cellfun(@(x) size(x,1), trialInds_byDir_byRew);
n_mincond = min(n_byDir_byRew(:));


predDirLabelsSS_byRew_byRep = cell(nrewards,nreps_mean);
decodeAccSS_byRew_byRep = nan(nrewards,nreps_mean);
fullStart = tic;

disp('Running normal subsampled CV for means')
for r = 1:nrewards
    tic
    curTrialInds_cells = trialInds_byDir_byRew(:,r);
    for i = 1:nreps_mean
        % Subsample, then get folds, labels, and data for the current rep
        repInds_cells = cellfun(@(x) x(randsample(length(x),n_mincond,false)), curTrialInds_cells,'uniformoutput',false);
        repInds_mat = [repInds_cells{:}];
        repFoldLabels_mat = repmat(crossvalind('kfold',n_mincond,nfolds),[1,ndirections]);
        repInds = repInds_mat(:);
        repFoldLabels = repFoldLabels_mat(:);
        repLabels = nan(n_mincond*ndirections,1);
        for d = 1:ndirections
            repLabels(((d-1)*n_mincond+1):d*n_mincond) = directions(d);
        end; clear d
        repData = binnedData(repInds,:);
        
        % Do normal CV decoding from here
        predDirLabelsSS_byRew_byRep{r,i} = quickerCVNBDecode(repData,repLabels,repFoldLabels);
        decodeAccSS_byRew_byRep(r,i) = 100*sum(predDirLabelsSS_byRew_byRep{r,i}==repLabels)./length(repLabels);
        if mod(i,500)==0,disp(['rep ' num2str(i) ', rew ' num2str(r)]); end
    end; clear i
    toc
end; clear r
toc(fullStart)

decodeAcc_byRew_mean = mean(decodeAccSS_byRew_byRep,2);

%% Next, do nested CV
MSE_byRew = nan(nrewards,1);

tic
disp('Running nested CV for error')
for r = 1:nrewards
    % Get relevant indices
    curInds = rewardLabels==rewards(r);
    curInput = binnedData(curInds,:);
    [~,curDirectionLabels] = ismember(directionLabels(curInds),directions); % reduces labels to 1:ndirections
    curDirCounts = histcounts(curDirectionLabels);
    n_mincond = min(curDirCounts);
    
    % Preallocate error lists
    a_list = nan(nreps_sem,nfolds); % (mean(inner)-mean(outer))^2 per fold/rep
    b_list = nan(nreps_sem,nfolds); % variance of outer per fold/rep
    es = nan(nreps_sem*(nfolds-1)*n_mincond*ndirections,1); % per rep, we accumulate (nfolds-1)*ndata here
    normalERRCV_byRep = nan(nreps_sem,1);        % normal CV'd error
    
    % Run it...
    for rep = 1:nreps_sem
        % pick indices for folds balanced-ly
        repInds = nan(ndirections,n_mincond);
        for d = 1:ndirections
            tempInds = find(curDirectionLabels==d);
            repInds(d,:) = tempInds(randsample(length(tempInds),n_mincond,false)); % pick n_mincond inds w/o replacement
        end; clear d
        
        % split into K folds
        repFoldLabels = crossvalind('kfold',n_mincond,nfolds);
        repOuterTestLabels = nan(numel(repInds),1);
        repOuterPredLabels = nan(size(repOuterTestLabels));
        for k = 1:nfolds
            % run crossval on non-k folds
            innerRepFoldInds = repInds(:,~ismember(repFoldLabels,k)); % anything not in fold K is fair game
            innerTestLabels = curDirectionLabels(innerRepFoldInds(:));
            innerFoldLabels = repmat(repFoldLabels(~ismember(repFoldLabels,k))',[ndirections 1]);
            innerPredLabels = quickerCVNBDecode(curInput(innerRepFoldInds(:),:),innerTestLabels,innerFoldLabels(:));
            innerErrors = innerTestLabels~=innerPredLabels; % this is e^{(in)}
            es(find(isnan(es),1)+(0:(length(innerErrors)-1))) = innerErrors;
                        
            % fit model with non-k folds and test on fold k. This is
            % basically just the result from normal CVing but for this
            % specific fold.
            trainRepInds = ~ismember(repFoldLabels,k);
            trainDataInds = flat(repInds(:,trainRepInds));
            trainData = curInput(trainDataInds,:);
            trainLabels = curDirectionLabels(trainDataInds);
            %             nbModel = fitcnb(trainData,trainLabels); % default for continuous values is gaussian NB; slow mode
            classMeans = cell(ndirections,1);
            classVars = cell(ndirections,1);
            for c = 1:ndirections
                classInds = trainLabels==c;
                classMeans{c} = mean(trainData(classInds,:));
                classVars{c} = var(trainData(classInds,:));
            end; clear c

            testRepInds = ismember(repFoldLabels,k);
            testDataInds = flat(repInds(:,testRepInds));
            testData = curInput(testDataInds,:);
            outerTestLabels = curDirectionLabels(testDataInds);
            LLs_cell = cellfun(@(x,y) -1/2*(sum(log(y)))-1/2*diag(((testData-x)*diag(1./y)*(testData-x)')), classMeans,classVars,'uniformoutput',false);
            LLs_mat = [LLs_cell{:}];
            [~,outerPredLabels] = max(LLs_mat,[],2);
            %             outerPredLabels = predict(nbModel,testData); %           slow mode
            outerErrors = outerTestLabels~=outerPredLabels; % this is e^{(out)}

            % Append normal cross-val fold result to the repOuter
            outerInds = find(isnan(repOuterTestLabels),1)+(0:(length(testDataInds)-1));
            repOuterTestLabels(outerInds) = outerTestLabels;
            repOuterPredLabels(outerInds) = outerPredLabels;
            
            % append all results
            a_list(rep,k) = (mean(innerErrors)-mean(outerErrors))^2;
            b_list(rep,k) = mean(outerErrors,1);  % normalize by N, not N-1
        end; clear k
        
        % Evaluate normal CVing results
        normalERRCV_byRep(rep) = mean(repOuterTestLabels~=repOuterPredLabels);
        if mod(rep,10)==0, disp(['Rew' num2str(r) ' rep' num2str(rep)]); end
    end; clear rep
    
    % Calculate stuff
    MSE_byRew(r) = mean(a_list(:))-var(b_list(:),1);
    r
end; clear r
toc

decodeAcc_byRew_sem = 100*sqrt((nfolds-1)/nfolds)*sqrt(MSE_byRew);

end

