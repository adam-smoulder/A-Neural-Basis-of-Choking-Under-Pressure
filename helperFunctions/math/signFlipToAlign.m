function [alignedB, signs] = signFlipToAlign(A,B)
% Takes matrix A that is nobs x ndims, and B which is ndims x nobs. We
% assume that A and B are the same (or nearly so), but the sign of the
% columns may be wrong. We find "signs" such that the dimensions of B match
% with A. The A matrix should be ndims x 1 as the reference to align to.
%
% Adam Smoulder, 9/23/21

A_signs = sign(A);
A_weights = A.^2; % higher vals should have bigger weight and be unsigned
% A_weights = sum(B.^2,2);

B_signs = sign(B);
simSign = abs(A_weights'*(A_signs-B_signs)); % should be 0 if they're the same sign
diffSign = abs(A_weights'*-(A_signs+B_signs)); % should be 0 if they're not the same sign
[~,flipInd] = min([simSign ; diffSign]); % if same sign, should be 1. If diff sign, should be 2.
signs = -2*flipInd+3; % If flipind is 2, this is -1. If flipind is 1, this is 1.
alignedB = B.*signs;
end

