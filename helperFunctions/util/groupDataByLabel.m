function [data_grouped,inds_grouped,n_grouped] = groupDataByLabel(data,labels)
% This function takes in some data where the first dimension is a number of
% trials (ntrials) and each trial has some label(s) associated with it. It
% outputs a cell matrix with the observations grouped by labels.
%
% Inputs:
% - data: [ntrials x __] matrix of data for each trial
% - labels: [ntrials x nlabels] matrix with the labels for each trial. Say
% we had a direction, reward size, and success label; one row may be [3 4
% 0] to indicate direction 3, reward 4, success 0.
%
% Outputs:
% - data_grouped: cell matrix with each dimension corresponding to one
% label (in the order of the columns of the labels input) of size equal to
% the number of unique values for each label. A given cell contains [nobs x
% __] data for each trial with that set of labels
% - inds_grouped: cell matrix with the indices corresponding to
% data_grouped
% - n_grouped: matrix with the length of each of the cells
%
% Adam Smoulder, 1/11/21 (edit 6/8/21 for adding n_grouped)

[ntrials,nlabels] = size(labels);
assert(ntrials==size(data,1))
assert(sum(isnan(sum(labels,2)))==0,'Nans in labels')
labels = double(labels); % in case it's boolean

% In a terrible design pattern, we use nans for preallocating for non-cell
% stuff, so if there are nans in the data, they'd get lost...so instead
% we're going to set these values to a specific magic number that we can
% re-set to nans later
if ~iscell(data)
    nanTempVal = 66666;
    data(isnan(data)) = nanTempVal;
end

% First, get all of the unique label values and convert the labels matrix
% to use indices instead of values
labelValues = cell(nlabels,1);
for a = 1:nlabels
    labelValues{a} = unique(labels(:,a));
    [~,labels(:,a)] = ismember(labels(:,a),labelValues{a});
end; clear a
nuniquelabels = cellfun(@(x) length(x), labelValues)';

% For each trial, store the values in the appropriate labels' cell.
% Regardless of how many labels/label-values there are, since we now
% have made the labels indexed (i.e., if there are 3 labels, they are
% always now [1,2,3]), we can just use:
% [label rows - 1] * [1 nlabel1 prod(nlabel1,nlabel2) . . .] + 1.
% Weird, but it helps circumvent differing numbers of labels!
% Further, since we don't know the dimensions of data beyond the first
% being ntrials, we'll just accumulate as we go by stacking rows. This
% isn't the most efficient way to do it with many trials, but is the
% easiest way to generalize to an unknown-sized data tensor
if length(nuniquelabels)~=1
    data_grouped = cell(nuniquelabels);
    inds_grouped = cell(nuniquelabels);
else
    data_grouped = cell(nuniquelabels,1);
    inds_grouped = cell(nuniquelabels,1);
end
labelMultiplier = [1 cumprod(nuniquelabels(1:end-1))]';
dataSize = size(data);

% % Method 1: Simpler but slower
% for i = 1:ntrials
%     cellIndex = (labels(i,:)-1)*labelMultiplier+1;
%     data_grouped{cellIndex}(end+1,:,:,:,:,:,:,:) = reshape(data(i,:),[1 dataSize(2:end)]);
%     inds_grouped{cellIndex}(end+1) = i;
% end; clear i


% Method 2: Preallocates for speed but is much uglier
for label = 1:numel(data_grouped)
    if ~iscell(data)
        data_grouped{label} = nan([ntrials dataSize(2:end)]);
    else
        data_grouped{label} = cell([ntrials dataSize(2:end)]);
    end
    inds_grouped{label} = nan([ntrials 1]);
end; clear label
for i = 1:ntrials
    cellIndex = (labels(i,:)-1)*labelMultiplier+1;
    data_grouped{cellIndex}(i,:,:,:,:,:,:,:) = reshape(data(i,:),[1 dataSize(2:end)]);
    inds_grouped{cellIndex}(i) = i;
end; clear i

if ~iscell(data)
    data_grouped = cellfun(@(x) x(~isnan(x(:,1)),:,:,:,:,:,:,:),data_grouped,'uniformoutput',false);
    data_grouped = cellfun(@(x) setValToNan(x,nanTempVal),data_grouped,'uniformoutput',false); % have to re-set those vals to nan
else
    data_grouped = cellfun(@(x) x(~cellfun(@(y) isempty(y),x),:,:,:,:,:,:),data_grouped,'uniformoutput',false);
end
inds_grouped = cellfun(@(x) x(~isnan(x(:,1)),:),inds_grouped,'uniformoutput',false);



n_grouped = cellfun(@(x) size(x,1), data_grouped);
end

function [input] = setValToNan(input,val)
input(input==val) = nan;
end










