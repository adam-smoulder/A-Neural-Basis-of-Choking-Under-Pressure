function [ellipsePts] = plotCovEllipse(mu,sigma,color,varargin)
%
% Plot a one-standard-deviation ellipse for a 2-dimensional Gaussian distribution
%
% INPUTS:
% mu - mean of the Gaussian (N x 1) 
% sigma - covariance matrix of the Gaussian (N x N)
% color - triple with color of contour
% optional input: p-value 
% 
% to be clear, p-value is "what % of the data should be covered by this
% ellipse (per a gaussian distribution with the input covariance)?"
%
% Code adapted from EM_GM.m by Patrick P. C. Tsui.

if nargin > 3
    pval = varargin{1};
else
%     % These are values for SD IFF you're in 1D
%     pval = 0.6827; % 1 standard codeviation 
%     pval = 0.95; % 2 standard codeviations
%     pval = 0.99; % 3 standard codeviations

    % In 2D
    pval = 0.3836; % In 2D, one standard dev. of each dim covers a total of 39% of the data
    % For a primer, check out this post 
    % http://johnthemathguy.blogspot.com/2017/11/statistics-of-multi-dimensional-data.html
    % I used to have a really good resource somewhere that explained this
    % well but can't find it...but if you sanity check do the following, it
    % should work:
    % temp = randn([10000,2]).*[1 4]; figure; plotCovEllipse(mean(temp)',cov(temp),[0 0 0])
    % So x-axis should go [-1 1] and y should go [-4 4], one std for each dir!
    %
    % alternatively, you can use the 1D interpretation of standard
    % deviation (I want an ellipse that covers ~69% of the data) and just
    % use that value instead. 
end

if size(mu,1) < size(mu,2), mu = mu'; end % make sure it's oriented correctly...

s = -2*log(1-pval);
[eigVecs, eigVals] = eig(s*sigma);
theta = linspace(0,2*pi,1000);
ellipsePts = (eigVecs*sqrt(eigVals))*[cos(theta); sin(theta)] + mu;
hold on
plot(ellipsePts(1,:),ellipsePts(2,:),'color',color,'linewidth',1.5)
% plot(mu(1),mu(2),'+','markeredgecolor',color,'markerfacecolor',color,'markersize',10,'linewidth',1)

end