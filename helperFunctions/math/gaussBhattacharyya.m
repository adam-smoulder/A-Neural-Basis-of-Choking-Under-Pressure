function [bhatDist, bhatCoeff] = gaussBhattacharyya(mean1, mean2, cov1, cov2)
% This function finds the Bhattacharyya distance between two multivariate
% distributions under the assumption that they're Gaussian distributed. 
%
% Inputs:
% - mean1 (and mean2): D x 1 vector representing the multidimensional mean of
% each distribution
% - cov1 (and cov2): D x D matrix of the covariance for the distribution
%
% Outputs:
% - bhatDist: Bhattacharyya distance between the distributions
% - bhatCoeff: the Bhattacharyya coefficient between the distributions
% (approximately measures overlap between samples)
%
% Adam Smoulder, 6/24/19

cov12 = (cov1+cov2)/2;
bhatDist = 1/8*(mean1-mean2)'/cov12*(mean1-mean2)+1/2*log(det(cov12)/(sqrt(det(cov1)*det(cov2)))); % straight from wikipedia...
bhatCoeff = exp(-bhatDist);

end

