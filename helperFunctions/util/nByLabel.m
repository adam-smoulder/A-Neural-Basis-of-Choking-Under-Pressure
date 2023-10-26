function [n_byLabel] = nByLabel(labels)
% This function finds the number of repetitions for each "condition" based
% on input labels. For instance, say each of 5 trials has a direction and
% reward label (in the 1st and 2nd column respectively) as follows:
% [3 2;
%  2 1;
%  1 1;
%  3 1;
%  3 2;]
% Since we have 3 unique direction labels and 2 unique reward labels, the
% output matrix would be [3x2]. We count each of the multiplexed
% conditions; hence the output would be:
% [1 0;  % where dir = 1 rew = 1, then where dir = 1 rew = 2
%  1 0;  % where dir = 2 rew = 1, then where dir = 2 rew = 2
%  2 1;] % where dir = 3 rew = 1, then where dir = 3 rew = 2
%
% Inputs:
% - labels: [ntrials x nlabels] matrix with the labels for each trial. Say
% we had a direction, reward size, and success label; one row may be [3 4
% 0] to indicate direction 3, reward 4, success 0.
%
% Outputs:
% - n_byLabel: matrix with each dimension corresponding to one label(in the 
% order of the columns of the labels input) with the number of times that
% condition showed up.
%
% NOTE labels are sorted lowest-to-highest.
%
% Adam Smoulder, 1/13/21

[ntrials,nlabels] = size(labels);

% First, get all of the unique label values and convert the labels matrix
% to use indices instead of values
labelValues = cell(nlabels,1);
for a = 1:nlabels
    labelValues{a} = unique(labels(:,a));
    [~,labels(:,a)] = ismember(labels(:,a),labelValues{a});
end; clear a
nuniquelabels = cellfun(@(x) length(x), labelValues)';

% For each trial, store the values in the appropriate labels' index.
% Regardless of how many labels/label-values there are, since we now
% have made the labels indexed (i.e., if there are 3 labels, they are
% always now [1,2,3]), we can just use:
% [label rows - 1] * [1 nlabel1 prod(nlabel1,nlabel2) . . .] + 1.
% Weird, but it helps circumvent differing numbers of labels!
n_byLabel = zeros(nuniquelabels);
labelMultiplier = [1 cumprod(nuniquelabels(1:end-1))]';
for i = 1:ntrials
    outputIndex = (labels(i,:)-1)*labelMultiplier+1;
    n_byLabel(outputIndex) = n_byLabel(outputIndex)+1;
end; clear i
end

