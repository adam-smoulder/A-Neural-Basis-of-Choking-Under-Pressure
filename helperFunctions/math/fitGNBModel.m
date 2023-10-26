function [modelMeans,modelVars,modelCovs] = fitGNBModel(data,labels)
% This function fits a Gaussian Naive Bayes model.
%
% Inputs:
% - data: nobs x ndims data input matrix
% - labels: nobs x 1 labels for nlabels different conditions
%
% Outputs:
% - modelMeans: nlabels x ndims means for each label
% - modelVars:  nlabels x ndims variances for each label
% - modelCovs:  nlabels x ndims x ndims covariance matrices for each label
% (just the diagonal of modelVars for each condition)
%
% Adam Smoulder, 9/27/21

groupedTrainData = groupDataByLabel(data,labels);
nlabels = length(unique(labels));
modelMeans = cellfun(@(x) mean(x), groupedTrainData, 'uniformoutput',false);
modelMeans = reshape([modelMeans{:}],[size(data,2) nlabels])';
modelVars = cellfun(@(x) var(x), groupedTrainData, 'uniformoutput',false);
modelVars = reshape([modelVars{:}],[size(data,2) nlabels])';
modelCovs = nan(nlabels,size(modelVars,2),size(modelVars,2));
for k = 1:nlabels
    modelCovs(k,:,:) = diag(modelVars(k,:));
end; clear k
end

