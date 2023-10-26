function [output] = sKLD(pmf1,pmf2)
% This function calculates the symmetric KL Divergence bewteen 2 pmfs. 
% We assume the PMFs are the same size. By "symmetric" I mean literally
% just taking the average of D(pmf1||pmf2) and D(pmf2||pmf1).
% 
% Adam Smoulder, 6/15/22
eps = 1E-10;
KLD1 = sum(pmf1.*log(pmf1./(pmf2+eps)+eps));
KLD2 = sum(pmf2.*log(pmf2./(pmf1+eps)+eps));
output = 0.5*(KLD1+KLD2);
end

