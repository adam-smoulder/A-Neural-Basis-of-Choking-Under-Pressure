function [shapeLabel,sigSigns,shapeNames] = evalRewardTuningShape3Rew(p,rewMeans)
% This function takes in the p-values and reward means associated with a
% Kruskal-Wallis test and attempts to classify the reward tuning in one of
% the following categories:
% - 1: Monotonic-up
% - 2: S-L up
% - 3: Monotonic-down
% - 4: S-L down
% - 5: U
% - 6: inv-U
% - 7: J-special up
% - 8: J-special down
% - 9: none
%
% To do this, we evaluate the significance and sign of changes in firing
% rate between the rewards and see how that shape matches one of the
% categories. FOR THIS FUNCTION we assume medium and large rewards have
% been combined! I'm just going to call it Large for this.
%
% Inputs:
% - p: [3 x 1] array of p-values, corresponding to the pairwise comparisons
% of each reward (SvL, SvJ, LvJ)
% - rewMeans: [3 x 1] array of mean FR for each reward
%
% Outputs:
% - shapeLabel: scalar, value from 1-9 corresponding to the different
% tuning shapes
% - sigSigns: [1 x 3] array of changes from S->L, S->J, and L->J, where -1 
% = downward, 0 = n.s., 1 = upward.
% - shapeNames: [9 x 1] cell array with the names associated with each of
% the 6 numerical labels
%
% Adam Smoulder, 12/7/22

% Get the signs associated with each comparison
diffMat = rewMeans-rewMeans';
signs = sign(diffMat([2 3 6])');
sigSigns = (p < 0.05).*signs;

% Now we're just going to painstakingly go over each possibility
shapeNames = {'mono-up','S-L-up','mono-down','S-L-down','U','inv-U','J-spec-up','J-spec-down','none'};
if sigSigns(1)==1 % S-L increase
    if sigSigns(3)==1 % L-J increase
        shapeLabel = 1; % mono-up
    elseif sigSigns(3)==-1 % L-J decrease
        shapeLabel = 6; % inv-U
    else % L-J no sig change
        if sigSigns(2)==1 && diffMat(3,1)>0 && diffMat(3,2) > 0 % greater than small and (ns) greater than L
            shapeLabel = 1; % mono-up
        else
            shapeLabel = 2; % S-L-up
        end
    end
elseif sigSigns(1)==-1 % S-L decrease
    if sigSigns(3)==1 % L-J increase
        shapeLabel = 5; % U
    elseif sigSigns(3)==-1 % L-J decrease
        shapeLabel = 3; % mono-down
    else % L-J no sig change
        if sigSigns(2)==-1 && diffMat(3,1)<0 && diffMat(3,2)<0 % less than small and (ns) less than L
            shapeLabel = 3; % mono-down
        else
            shapeLabel = 4; % S-L-down
        end
    end
else % S-L no sig change
    if sigSigns(3)==1 % L-J increase
        shapeLabel = 7; % J-spec-up
    elseif sigSigns(3)==-1 % L-J decrease
        shapeLabel = 8; % J-spec-down
    else % L-J no sig change
        if sigSigns(2)==1 % S-J increase
            shapeLabel = 7; % J-spec-up
        elseif sigSigns(2)==-1 % S-J decrease
            shapeLabel = 8; % J-spec-down
        else % no sig changes at all
            shapeLabel = 9; % none
        end
    end
end
end

