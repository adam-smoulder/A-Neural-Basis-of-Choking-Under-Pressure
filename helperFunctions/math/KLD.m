function [output] = KLD(pmf1,pmf2)
% This function calculates the KL Divergence bewteen 2 pmfs. 
% We assume the PMFs are the same size. 
% 
% Adam Smoulder, 6/28/22
eps = 1E-12;
output = sum(pmf1.*log(pmf1./(pmf2+eps)+eps));
end