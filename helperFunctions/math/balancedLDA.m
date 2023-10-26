% INSERT DESCRIPTION HERE
%
function [w, objVal, eigVls, classMeans, sW] = balancedLDA(input, classLabels, groupLabels)

% Group data by condition
input_byClass_byGroup = groupDataByLabel(input,[classLabels groupLabels]);
n_byLabel_byGroup = cellfun(@(x) size(x,1), input_byClass_byGroup);
[nclasses,ngroups] = size(n_byLabel_byGroup);
ndims = size(input_byClass_byGroup{1,1},2);

% Get the means for each class x group (unbiased), then take the mean
% across groups, and calculate scatter-between (sB, also sometimes called
% signal covariance) from these class means. We do this instead of getting
% sB for each group as sB is biased with low sample counts
inputMeans_byClass_byGroup = cellfun(@(x) mean(x)', input_byClass_byGroup,'uniformoutput',false); % each cell is now ndims x 1
classMeans = nan(nclasses,ndims);
for c = 1:nclasses
    classMeans(c,:) = mean([inputMeans_byClass_byGroup{c,:}]');
end; clear c
sB = cov(classMeans,1); % normalize by 1/N for this, not 1/(N-1)

% Next, get the covariance matrix for each class x group, then average all
% of these equally. We calculate covariance in an unbiased manner, so the
% overall scatter-within matrix that averages over these (sW, also
% sometimes called noise covariance) should treat them all equally.
inputCovs_byClass_byGroup = cellfun(@(x) permute(cov(x),[1 3 2]), input_byClass_byGroup,'uniformoutput',false); % permute adds a singleton dimension to the middle for easy concatenation
sW = squeeze(mean([inputCovs_byClass_byGroup{:}],2));


% Finally, solve the generalized eigenvalue problem for sB/sW; this is core
% computation of LDA. The only difference from usual is that we got sB and
% sW treating all groups as equally contributing.
[eigVecs, eigVls] = eig(sB, sW); % these are in reverse order for some reason... let's fix that
[eigVls, order] = sort(diag(eigVls),'descend');
eigVecs = eigVecs(:,order);

% find the rank (# labels - 1 or #dims) and extract w accordingly
m = min(nclasses, ndims);
w = eigVecs(:, 1:(m-1));
objVal = det(w'*sB*w)/det(w'*sW*w);
eigVls = eigVls(1:(m-1));

end

