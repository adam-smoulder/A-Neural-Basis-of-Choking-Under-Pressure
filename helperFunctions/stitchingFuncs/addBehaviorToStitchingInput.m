function [stitchingInput] = addBehaviorToStitchingInput(stitchingInput,behavior)
% This function takes the stitchingInput structure and a post-analysis
% behavioral structure (which has reactionTime, goCueTime, etc) and
% combines the two into an output stitchingInput structure, where the
% behavioral fields have been individually inserted into the
% stitchingInput and are stretched from a trial-to-trial basis to a time
% bin to time bin label.
%
% Adam Smoulder, 5/28/19

% First, get the day*trial labels from masterInfo and from behavior
masterInfo = stitchingInput.masterInfo;
dayTrial_fromMaster = masterInfo(:,1)*1000+masterInfo(:,2);
dayTrial_fromBehavior = behavior.day*1000+behavior.trial;

% Now, stretch behavior's data based on the same day-trial labels
for i = 1:length(dayTrial_fromMaster) % we loop through it this way to appropriately place nans where data is missing
    curDayTrial = dayTrial_fromMaster(i);
    behInd = dayTrial_fromBehavior==curDayTrial;
    if sum(behInd)==1 % normal case, just one index matches
        vals = [behavior.endReachTime(behInd) ; behavior.goCueTime(behInd) ; ...
            behavior.initialSpeed(behInd) ; behavior.initialVelocity(behInd) ;  ...
            behavior.peakSpeed(behInd) ; behavior.reactionTime(behInd)]; % we use "vals" so I don't have to redefine everything for the conditional here
    elseif sum(behInd)>1 % odd case, we have multiple RTs for one trial?
        disp('ADAM PROBLEM PROBLEM'); pause
    else % we have RT and stuff but no data in master info for it; skip it
        vals = nan(6,1);
    end
    
    % and assign it!
    stitchingInput.endReachTime(i) = vals(1);
    stitchingInput.goCueTime(i) = vals(2);
    stitchingInput.initialSpeed(i) = vals(3);
    stitchingInput.initialVelocity(i) = vals(4);
    stitchingInput.peakSpeed(i) = vals(5);
    stitchingInput.reactionTime(i) = vals(6);
    
    % redundant, but put behavior in as well
    stitchingInput.behaviorByTrial = behavior;
    
    if mod(i,10000)==0, disp(['Completed ' num2str(i) ' time bins of assignment']); end
end; clear i

end

