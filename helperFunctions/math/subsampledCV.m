function [estMeans,estSEMs] = subsampledCV(X,Y,Z,method,K,R,varargin)
% This function performs subsampled cross-validation to ensure that
% mean decoding accuracies are not biased due to having greater amounts of
% training data.
%
% Inputs:
% - X: [N x D] data, where N is number of observations and D is number of
% dimensions. 
% - Y: [N x 1] labels that we are decoding from the data
% - Z: [N x 1] labels that we wish to balance trial counts over
% - method: string, the classification/regression method to be used.
% Options include 'GNB','LDA'
% - K: scalar, number of folds to use (must be >= 2)
% - R: scalar, number of repetitions to use.
%
% Optional inputs:
% - LOGGING: boolean; if true, will show update every 500 reps (default:
% true)
%
% Outputs:
% - estMean: [Cz x 1], the mean decoding accuracy for each of the Cz unique
% classes of Z.
% - estSEMs: [Cz x 1], the average of SEMs over folds for each class of
% Z. This is often viewed as a conservative metric for error (but with
% subsampling...who knows? I wouldn't use this. Use nestedCV.m within each
% of the Cz classes instead).
%
% As an example, I am using this because I want to decode the upcoming
% reach direction (Y, 1:4 or 1:8) from my neural data (X) in different
% reward conditions (Z, 1:4). In particular, one of the reward conditions 
% (4) has many less samples than the others, about 1/5 as many. So, I want
% to subsample my Z classes 1:3 down many times and get cross-validated
% decoding accuracy for reach direction such that they can be compared to
% results from Z = 4.
%
% Adam Smoulder, 10/4/21

if ~isempty(varargin)
    LOGGING = varargin{1};
else
    LOGGING = true;
end

% If KNN, get K for that
if contains(method,'NN') % KNN; get the K
    KNN_K = str2double(method(1:strfind(method,'NN')-1));
    method = 'KNN';
end

% Get overall size of data and group data
[~,inds_byLabel] = groupDataByLabel(X,[Y Z]);
[Cy,Cz] = size(inds_byLabel); % number of classes in Z
n_byLabel = cellfun(@(x) length(x), inds_byLabel);
n_mincond = min(n_byLabel(:)); % number of points to subsample to

% In clase the classes of Y aren't 1:Cy, replace them
Yvals = unique(Y);
[~,Y] = ismember(Y,Yvals);
Yvals = (1:Cy)';

% We will scramble our indices each time, so we only have to generate fold
% labels once
foldLabels = repmat(crossvalind('kfold',n_mincond,K),[Cy 1]);

% Run through each repetition for each class of Z
decodeAcc_byZ_byRep_byFold = nan(Cz,R,K);
for cz = 1:Cz % for each class of Z
    for r = 1:R % for each repetition
        % Scramble our indices within each class and subsample this rep's data
        repPerm = cellfun(@(x) randperm(size(x,1))', squeeze(inds_byLabel(:,cz)),'UniformOutput',false);
        repInds_cell = cellfun(@(x,y) x(y(1:n_mincond),:), squeeze(inds_byLabel(:,cz)), repPerm, 'UniformOutput', false);
        repInds = cellArrayToVector(repInds_cell);
        
        % Now do cross-validation
        for k = 1:K
            % Get train/test inds
            repTrainInds = foldLabels~=k; % indices in our subsample
            trainInds = repInds(repTrainInds); % indices in the full data
            repTestInds = foldLabels==k;
            testInds = repInds(repTestInds);
            testY = Y(testInds);
            
            % Train and test
            switch method
                case 'GNB'
                    [modelMeans,~,modelDiagCovs] = fitGNBModel(X(trainInds,:), Y(trainInds)); % train
                    predY = gaussMLE(X(testInds,:), modelMeans, modelDiagCovs);
                case 'LDA'
                    [~, ~, ~, ~,modelMeans,modelNoiseCov] = linearDiscAnalysis(X(trainInds,:), Y(trainInds),'DO_CV',false);
                    predY = gaussMLE(X(testInds,:), modelMeans, repmat(permute(modelNoiseCov,[3 1 2]),[Cy 1 1]));
                case 'KNN'
                    predY = nan(length(testInds),1);
                    for i = 1:length(testInds)
                        distVals = sum((X(trainInds,:)-X(testInds(i),:)).^2,2);
                        [~,order] = sort(distVals,'ascend');
                        closeY = Y(trainInds(order(1:KNN_K)));
                        count_byClass = histcounts(closeY,(0:Cy)+0.5);
                        [maxVal,predYVal] = max(count_byClass);
                        if sum(count_byClass==maxVal) > 1 % multiple classes with same; pick the one with closer values
                            validY = find(count_byClass==maxVal);
                            meanDist = inf;
                            for y = validY
                                curMeanDist = sum((closeY==y).*distVals(order(1:KNN_K)));
                                if curMeanDist < meanDist % If this is the closest yet, keep it
                                    predYVal = y;
                                    meanDist = curMeanDist;
                                end
                            end; clear y
                        end
                        predY(i) = predYVal;
                    end; clear i
                otherwise
                    error('Method not recognized')
            end
            decodeAcc_byZ_byRep_byFold(cz,r,k) = 100*mean(Yvals(predY)==testY); % we convert to percentage here
        end; clear k
        
        if mod(r,500)==0 && LOGGING, disp(['Completed rep ' num2str(r) ' for group ' num2str(cz)]); end
    end; clear r
end; clear cz

% Take the average and SE over all reps/folds
estMeans  = mean(decodeAcc_byZ_byRep_byFold(:,:),2);
estSEMs = mean(std(decodeAcc_byZ_byRep_byFold,[],3)/sqrt(K),2);

% Scale the SE by the actual numbers of trials
estSEMs = estSEMs.*sqrt(n_mincond*Cy./sum(n_byLabel)');





end

