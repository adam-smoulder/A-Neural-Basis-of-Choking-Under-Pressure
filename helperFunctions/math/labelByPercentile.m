function [prctileLabels,binEdges] = labelByPercentile(data,nclasses)
% This function takes a vector of data and gives class labels (1,2,3...)
% based on what percentile the data is in. 
%
% Inputs:
% -data: [Nx1] vector of data, assumed to be numbers
% -nclasses: scalar of how many classes you want to divide the data into
%
% Outputs:
% -prctileLabels: [Nx1] vector of labels
% -binEdges: [nclasses+1 x 1] vector of the edges used for each bin
%
% E.g.
% data = [1 3 2 5 6 2 4];
% nclasses = 3;
% labelByPercentile(data,nclasses)
% Result: [1 2 1 3 3 1 2]
%
% Adam Smoulder, 6/23/20

% Remove nans
data(isnan(data),:,:) = [];

prcStepToUse = 100/nclasses;
prcsToUse = prcStepToUse:prcStepToUse:100;
prctileLabels = ones(length(data),1);
binEdges = nan(nclasses+1,1);
binEdges(1) = min(data);
for i = 2:nclasses
    if i == nclasses
        prctileLabels(data >= prctile(data,prcsToUse(end-1))) = nclasses;
    else
        prctileLabels(data >= prctile(data,prcsToUse(i-1)) & data <  prctile(data,prcsToUse(i))) = i;
    end
    binEdges(i) = prctile(data,prcsToUse(i-1));
end; clear i
binEdges(end) = max(data);

end

