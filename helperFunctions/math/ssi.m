function [SSI,isp,binEdges,priorDist,likeDist,dataDist,postDist] = ssi(data,labels,varargin)
% This function calculates the stimulus specific information (SSI), a
% metric that represents how much of the information in each response of
% each predictor is accounted for by each unique stimulus label.
%
% Mathematically, we get this by first getting the specific information:
% 
% isp(r) = H(S) - H(S|R=r)
%
% where r = response, s = stimulus. This represents the reduction in
% uncertainty about the stimulus given a certain response. H is the shannon
% entropy (assumes the function "shannonEntropy.m" is in the filepath)
%
% From there, we then get the SSI by taking the conditional expectation
% of the isp for the given stimulus:
%
% SSI(s) = E[isp | S = s] = sum over i in all response bins of ( P(R=ri|S=s)*isp(ri))
%
% Inputs:
% - data: N x D matrix of data, where there are N observations and D
% dimensions. For SSI, each of these dimensions is looked at individually.
% This could be neurons, factors, whatever.
% - labels: N x 1 vector of stimulus labels. Here, we assume the stimulus
% is NOT continuous and there are K unique stimuli
%
% Optional Inputs (pass through as name value pairs)
% - 'nbins': integer, number of bins to use for data histograms. Uses the
% same number for all dimensions of data because I'm lazy and don't want to
% go into cell-land...just loop over this for individual dimensions if you
% want to use different numbers of bins. Default is 30
% - 'binEdges': D x nbins+1 matrix, edges for the bins to use for the data.
% Default is determined from the nbins and the range of the data.
%
% Outputs:
% - SSI: D x K matrix with each dimension's SSI for each stimulus
% - isp: D x nbins matrix with the specific information of each of the
% binned responses for each dimension
% - binEdges: D x nbins+1 matrix with the edges of the bins for each
% histogram that was made for the data's distributions
% - priorDist: 1 x K vector of prior probabilities for each stimulus
% - likeDist: D x K x nbins matrix with the likelihood distributions; this
% is the probability of the response given a certain stimulus for each dim
% - dataDist: D x nbins matrix with the data distributions for each
% dimension across the responses binned by binEdges
% - postDist: D x K x nbins matrix with the posterior distributions; this
% is the probability of a stimulus given a certain response for each dim
% 
% Adam Smoulder, 10/14/19 (edit 10/15/19)

% set our defaults
trueBinEdges=[]; % empty unless overwritten
nbins = 30;
for i = 1:length(varargin)
    if strcmp(varargin{i},'nbins')
        nbins = varargin{i+1};
    end
    
    if strcmp(varargin{i},'binEdges')
        trueBinEdges = varargin{i+1};
        nbins = length(trueBinEdges)-1;
    end
end

% get some quick stuff we'll need
[N,D] = size(data);
assert(N==length(labels));
stimVals = unique(labels);
K = length(stimVals);

% First, let's get the prior distribution
priorCounts = histcounts(labels,[stimVals(1)-1E-10 ; stimVals+1E-10]); % histcounts uses bin centers, so we offset just a smidget
priorDist = priorCounts/sum(priorCounts); % Should be K x 1

% Next, we are going to make our likelihood and posterior distributions. To
% do this, we first have to bin the data in a consistent manner, spread
% enough to have a good pmf but tight enough to get sufficient observations
% for each bin. Set nbins at the top to change this.
binEdges = nan(D,nbins+1);
likeDist = nan(D,K,nbins);
dataDist = nan(D,nbins);
for d = 1:D
    if isempty(trueBinEdges)
        curBinEdges = linspace(min(data(:,d))-1E-5, max(data(:,d))+1E-5,nbins+1);
    else
        curBinEdges = trueBinEdges(d,:);
    end
    binEdges(d,:) = curBinEdges;
    for k = 1:K % make data|stimulus histograms 
        curData = sort(data(labels==stimVals(k),d),'ascend');
        likeCounts = histcounts(curData,curBinEdges);
        likeDist(d,k,:) = likeCounts/sum(likeCounts);
    end; clear k
    
    % We'll need the data distribution for each dim to get the posterior
    dataCounts = histcounts(data(:,d),curBinEdges);
    dataDist(d,:) = dataCounts/sum(dataCounts);
end; clear d

% The posterior is likelihood*prior/data
postDist = likeDist.*priorDist./permute(dataDist,[1 3 2]); 
for d = 1:D % 0/0 = nans; these should be replaced with the prior
    for i = 1:nbins
        if isnan(sum(postDist(d,:,i)))
            postDist(d,:,i) = priorDist;
        end
    end; clear i
end; clear d
assert(sum(sum(abs(sum(postDist,2)-1)<1E-10)) == nbins*D) % assert they all sum to 1

% isp is the entropy of the prior minus the posterior ; H(S)-H(S|X)
isp = nan(D,nbins);
for d = 1:D
    for i = 1:nbins
        isp(d,i) = shannonEntropy(priorDist)-shannonEntropy(postDist(d,:,i));
    end; clear i
end; clear d

% VERY finally, SSI is the conditional expectation of isp across responses
SSI = nan(D,K);
for k = 1:K
   SSI(:,k) = sum(squeeze(likeDist(:,k,:)).*squeeze(isp),2);
%      SSI(:,k) = sum(squeeze(likeDist(:,k,:).*permute(isp,[1 3 2]))); % weird permute is to avoid issues with squeeze on likeDist for D = 1 case
%    SSI(:,k) = sum(squeeze(likeDist(:,k,:).*permute(isp,[1 3 2])),2); % weird permute is to avoid issues with squeeze on likeDist for D = 1 case

end; clear k

end

