function [shapeLabel,changes,shapeNames] = evalRewardTuningShape(p,rewMeans)
% This function takes in the p-values and reward means associated with a
% Kruskal-Wallis test and attempts to classify the reward tuning in one of
% the following categories:
% - 1: Monotonic-up
% - 2: Monotonic-down
% - 3: U
% - 4: inverted-U
% - 5: other
% - 6: none
% To do this, we evaluate the significance and sign of changes in firing
% rate between the rewards and see how that shape matches one of the
% categories.
%
% Inputs:
% - p: [6 x 1] array of p-values, corresponding to the pairwise comparisons
% of each reward (SvM, SvL, SvJ, MvL, MvJ, LvJ)
% - rewMeans: [4 x 1] array of mean FR for each reward
%
% Outputs:
% - shapeLabel: scalar, value from 1-6 corresponding to the different
% tuning shapes
% - changes: [1 x 3] array of changes from S->M, M->L, and L->J, where -1 =
% downward, 0 = n.s., 1 = upward. In certain cases, these are changed by
% further comparisons (e.g. if S->M is 0 and M->L is 0, L will be compared
% to S, and if sig, will change the M->L value to the S->L change). These
% give an idea of the shape of the tuning curve
% - shapeNames: [6 x 1] cell array with the names associated with each of
% the 6 numerical labels
%
% Assumes 4 rewards are used.
%
% Adam Smoulder, 8/3/22

% Get the changes for adjacent rewards
changes = ((p([1 4 6])<0.05).*sign(diff(rewMeans)))';

% If 2x n.s. in a row, compare to previous value(s)
doubleSame = strfind(changes,[0 0]);
for ind = 1:length(doubleSame)
    if doubleSame(ind)==1 % SS? = look at large vs small
        changes(2) = (p(2)<0.05)*sign(diff(rewMeans([1 3])));
    elseif doubleSame(ind)==2 && changes(2)==0  % ?SS = look at J vs M
        changes(3) = (p(5)<0.05)*sign(diff(rewMeans([2 4])));
        if changes(3)==0 && changes(1)==0 % if still unresolved and no sig change form S, try S v J
            changes(3) = (p(3)<0.05)*sign(diff(rewMeans([1 4]))); % This is the really specific case of a monotonic up/down, where the only sig is S vs. J
        end
    end
end; clear ind

% Based on the changes, select an output label. We first have to declare
% the case for each:
shapeNames = {'mono-up','mono-down','U','inv-U','other','none'};
labelClasses = {... % each cell has cells containing all valid changes for that label
    {[1 1 1],[1 1 0],[1 0 1],[1 0 0],[0 1 1],[0 1 0],[0 0 1]};
    {[0 0 -1],[0 -1 0],[0 -1 -1],[-1 0 0],[-1 0 -1],[-1 -1 0],[-1 -1 -1]};
    {[0 -1 1],[-1 1 1],[-1 1 0],[-1 0 1],[-1 -1 1]};
    {[1 1 -1],[1 0 -1],[1 -1 0],[1 -1 -1],[0 1 -1]};
    {[-1 1 -1],[1 -1 1]};
    {[0 0 0]};
    };
% Now see which the current changes are in
shapeLabel = find(cellfun(@(y) sum(cellfun(@(x) sum(x==changes)==3,y)), labelClasses));
end

