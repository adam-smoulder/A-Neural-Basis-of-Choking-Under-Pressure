function [LL] = gaussLL(x,mu,sigma,varargin)
% For datapoint(s) x, gets the log-likelihood for multivariate gaussian
% distribution based on parameters mu and sigma.
%
% Inputs:
% - x: [npts x ndims] data to evaluate
% - mu: [1 x ndims] mean of gaussian distribution
% - sigma: [ndims x ndims] covariance of gaussian distribution
%
% Optional Inputs:
% - invSigma: [ndims x ndims] inverse of covariance matrix. Having this
% provided massively speeds up things for getting LL trial by trial
%
% Outputs:
% - LL: [npts x 1] log-likelihood of each datapoint
%
% Validated using the following commands (should produce same results):
% mu = zeros(1,10);
% sigma = eye(10);
% x = mvnrnd(mu,sigma,20);
% [gaussLL(x,mu,sigma) log(mvnpdf(x,mu,sigma))]
%
% Note - this equivalence may not work for higher dimensionality, as then
% the output of mvnpdf may get below machine precision. Hence why I made
% this...
%
% Adam Smoulder, 9/20/21

if nargin > 3
    invSigma = varargin{1};
else
    invSigma = inv(sigma);
end

[N,D] = size(x); % ndims

term1 = -D/2*log(2*pi);
term2 = -1/2*logdet(sigma);
if N < 5000 % faster to use linear algebra and use diagonal (no idea what the actual right number is, diag is faster at least through N = 1000 though)
    term3 = -1/2*diag((x-mu)*invSigma*(x-mu)');
else % faster to just do each individually
    term3 = nan(N,1);
    for n = 1:N
        term3(n) = -1/2*(x(n,:)-mu)*invSigma*(x(n,:)-mu)';
    end; clear n
end

LL = term1+term2+term3;

end

