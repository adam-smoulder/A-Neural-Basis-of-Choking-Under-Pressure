function [corrval] = nancorr2var(data)
% Finds the correlation matrix between variables while ignoring nans. Not
% sure why this doesn't exist...the present implementation is certainly not
% optimized, but it works and doesn't take too long (a few seconds for a
% 10000 x 250 dataset for me)
%
% THIS is just a special version of my normal nancorr, but just getting it
% bewteen two variables to spit out a scalar
%
% If two variables have no overlap, their correlation will be returned as
% nan.
%
% Inputs:
% - data: N x 2 matrix for N observations of 2 variables, where unobserved
% instances are filled in with nans
%
% Outputs:
% - corrVal: scalar, the correlation between the two series
%
% Adam Smoulder, 1/26/21
nanInds = isnan(sum(data,2));
if sum(~nanInds) <=2
    corrval = nan;
    return
end
data(nanInds,:) = [];
[N,K] = size(data);
indivVar = nanvar(data)'; % individual variances excluding nans

% We make a covariance matrix, excluding nans
covmat = zeros(K,K); % we use zeros to add the diagonal for indivVar
covmat = covmat+diag(indivVar);
for i = 1:(K-1)
    for j = (i+1):K
        curData = data(:,[i j]);                    % get current data
        curData(isnan(sum(curData,2)),:) = [];      % remove observations where either variable is nan
        if isempty(curData)                         % if there was no overlap b/w the variables
            covmat(i,j) = nan; covmat(j,i) = nan;   % the cov is nan
        else                                        % otherwise, we use the covariance of the observations remaining
            curCov = cov(curData);                  % though we don't need the full pair matrix, so just get the covariance value!
            covmat(i,j) = curCov(1,2); covmat(j,i) = curCov(1,2);
        end
    end; clear j
    % i
end; clear i

% finally, to get correlation, we simply divide the covariance by the
% square root of the individual variances of each pair
denominatorMatrix = sqrt(indivVar*indivVar');              % should be K x K with the product of each variables's stdev
corrmat = covmat./denominatorMatrix;
corrval = corrmat(1,2);

if abs(corrval) > 1
    disp('er wtf')
end
end

