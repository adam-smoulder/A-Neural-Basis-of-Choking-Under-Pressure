function [Q_t,Q_all] = tanglingCalc(trajMat,dt,varargin)
% Calculates tangling per Russo et al. 2018. 
%
% Inputs:
% - trajMat: [N x T] matrix of N dimensions (e.g. neurons/latents) as rows,
% each with T observations. Should be a single (typically mean) trajectory
% across the dimensions.
% - dt: scalar that is the amount of time between observations (s).
%
% Optional Inputs:
% - eps: scalar value to use in denominator to ensure noninfinite
% calculations; this also effectively becomes the minimum value of
% tangling, so it should be consistent across different conditions.
% 
% Outputs:
% - Q_t: [1 x T-2] vector with the max tangling for each time-point along
% the trajectory. 2 timepoints are lost to taking and realigning the
% derivative.
% - Q_all: [T-2 x T-2] full tangling matrix; the "un-maxed" version of Q_t
%
% Adam Smoulder, 6/30/20
[~,T] = size(trajMat);
time = 0:dt:(T-1)*dt;

% Calculate derivative & realign
d_trajMat_unaligned = diff(trajMat,[],2)./dt;
d_trajMat = twoPointInterpolation(time,...  % The deriv is most accurate to the midpoint of our time steps (e.g. we do the difference
   time-dt/2, time+dt/2,...             % b/w T = 1 and 2, so the deriv is most accurate to T = 1.5), so we linearly interpolate between
   d_trajMat_unaligned(:,1:end-1),d_trajMat_unaligned(:,2:end)); % deriv values to align back to the original T
trajMat = trajMat(:,2:end-1); % shave off start/end to match derivative

% Calculate tangling (note it's symmetric, I'm just being lazy/inefficient)
Q_num_all = nan(T-2,T-2); % (t,t') (not that it matters) - numerator of Q_all
Q_denom_all = nan(T-2,T-2); % (t,t') - denominator of Q_all
if ~isempty(varargin)
    eps = varargin{1};
else
    eps = 0.1*mean(sum(trajMat.^2)); % epsilon, calculated as done in the paper
end
for t = 1:(T-2)
%     Q_all(t,:) = sum((d_trajMat(:,t)-d_trajMat).^2)./(sum((trajMat(:,t)-trajMat).^2)+eps);
    Q_num_all(t,:) = sum((d_trajMat(:,t)-d_trajMat).^2);
    Q_denom_all(t,:) = sum((trajMat(:,t)-trajMat).^2)+eps;
end; clear t 
Q_all = Q_num_all./Q_denom_all;
Q_t = max(Q_all);
end

