function [pvals] = binomialProportionTest2(y,n)
% This function performs a binomial proportion test between the proportions
% in the vector y with counts n. This is the same test as
% binomialProportionTest.m but generalized to just take the proportion
% (i.e., the fraction of times the event happened) and n. This is useful in
% cases where the proportion isn't strictly integer (e.g., the result of
% cross-validated decoding accuracy).
%
% Adam Smoulder, 8/31/22

assert(length(y)==length(n),'Inputs are different size')
assert(((size(y,1)==1) || size(y,2)==1) && length(size(y))==2,'y is not a vector')
if any(y > 1) && sum(y <= 100)==length(y) % seems to be %
    disp('Binomial proportion test, converting percentages to fractions')
    y = y/100;
elseif any(y > 1) && all(y < n) % liekly counts; convert to fractions
    disp('Binomial proportion test, converting counts to fractions')
    y = y./n;
end

% Get the z-statistic and p-val for each pair
phat = (n.*y+(n.*y)')./(n+n'); % pooled success rates
z = (y-y')./sqrt(phat.*(1-phat).*(1./n+1./n'));
pvals = min(normcdf(z),normcdf(-z))*2; % taking the min and *2 bc two-tailed

end

