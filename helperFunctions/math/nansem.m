function [Y] = nansem(X,DIM)
% Calculates standard error of the mean, ignoring nan values.
%
% Adam Smoulder, 7/14/20 (edited 3/27/23)

Y = nanstd(X,[],DIM)./sqrt(sum(~isnan(X),DIM));
if isempty(Y), Y = nan; end

end

