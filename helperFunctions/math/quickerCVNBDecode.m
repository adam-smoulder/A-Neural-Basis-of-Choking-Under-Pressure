function [predLabels] = quickerCVNBDecode(inputData,inputLabels,foldLabels)

% Normal slow matlab implementation (variable names are diff...)
%         predDirLabelsSS_byRew_byRep{r,i} = nan(length(repFoldLabels),1);
%         for f = 1:nfolds
%             trainInds = repFoldLabels~=f;
%             trainData = repData(trainInds,:);
%             trainLabels = repLabels(trainInds);
%             nbModel = fitcnb(trainData,trainLabels);
%             
%             testInds = repFoldLabels==f;
%             testData = repData(testInds,:);
%             testLabels = repLabels(testInds);
%             predDirLabelsSS_byRew_byRep{r,i}(testInds) = predict(nbModel,testData);
%         end; clear f

% Might be faster than MATLAB's implementation...
origFolds = unique(foldLabels);
nfolds = length(origFolds);
predLabels = nan(size(foldLabels));
origClasses = unique(inputLabels);
nclasses = length(origClasses);
[~,inputLabels] = ismember(inputLabels,origClasses); % convert labels to 1:nclasses instead of whatever they are
for f = 1:nfolds
    % Get means and variances for GNB model; assume uniform prior
    trainInds = foldLabels~=origFolds(f);
    trainData = inputData(trainInds,:);
    trainLabels = inputLabels(trainInds);
    classMeans = cell(nclasses,1);
    classVars = cell(nclasses,1);
    for k = 1:nclasses
        curInds = trainLabels==k;
        classMeans{k} = mean(trainData(curInds,:));
        classVars{k} = var(trainData(curInds,:));
    end; clear k
    
    % Test; assumes number of dimensions is same across obs
    testInds = find(foldLabels==origFolds(f));
    testData = inputData(testInds,:);
    
    % Option 1: do it in bulk; should be faster for normal/low
    % dimensionality (i.e. 7x faster than the trial-by-trial option for D =
    % 20)
    LLs_cell = cellfun(@(x,y) -1/2*(sum(log(y)))-1/2*diag(((testData-x)*diag(1./y)*(testData-x)')), classMeans,classVars,'uniformoutput',false);
    LLs_mat = [LLs_cell{:}];
    
%     % Option 2: trial by trial
%     LLs_mat = nan(length(testData),nclasses);
%     for i = 1:length(testData)
%         LLs_mat(i,:) = cellfun(@(x,y) -1/2*(sum(log(y)))-1/2*((testData(i,:)-x)*diag(1./y)*(testData(i,:)-x)'), classMeans,classVars);
%     end; clear i
    
    [~,predLabels(testInds)] = max(LLs_mat,[],2);
    
end; clear f
predLabels = origClasses(predLabels); % convert back to the original class convention

end

