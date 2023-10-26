function [trialStatusLabels] = combineSimilarTrialStatuses(trialStatusLabels,version)
% This function combines statuses of similar nature into the same label.
% For instance, we have 3 classes of overshoots with different labels: -22
% (miss), -31 (scuff), and -32 (blow-through). These can all be set to the
% same value (like -22).
%
% We also have certain statuses that we want to ignore or exclude often 
% that can all be set to the same label. For instance, early returns (-33)
% and wild (-20) almost never happen, so we should just probably remove
% these trials, along with quitouts (0), so we can set them all to 0.
%
% Inputs:
% - trialStatusLabels: [ntrials x 1] vector with the status label for each
% trial.
% - version: string with what setting we want to use for combining.
%
% Outputs:
% - trialStatusLabels: same as input, but with certain labels combined.
%
% The full listing of trial statuses and their meaning as of 1/13/21 is at
% the bottom of this function.
%
% Adam Smoulder, 1/13/21 (edit 12/8/22)

if strcmp(version,'2020_behaviorPaper')
    oldStatuses = {[0 -13 -14 -20 -21 -33]; ...
                   [-22 -31 -32]; ...
                   [-23 -24 -25]; ...
                   [-34 -35]};
    newStatuses = [0; -22; -23; -34];
end

% For each entry in oldStatuses, set alllll of those labels to the
% corresponding index in newStatuses
for i = 1:length(oldStatuses)
    curInds = ismember(trialStatusLabels,oldStatuses{i});
    trialStatusLabels(curInds) = newStatuses(i);
end; clear i



% Status label meanings:
%   0 = Quitout: an active reach away from the direction of the reach target
% -11 = False Start: a reach attempt was made towards the targetbefore the go cue had been presented
% -12 = Delay Drift: hand lazily drifted out of the center, sometimes towards the target ("cheat") others not
% -20 = Wild: reach is too far off the course to call it either an overshoot or undershoot
% -21 = No Attempt: no attempt (or a very lazy one) was made, indicated by a very low or super late peak speed
% -22 = Overshoot (type 1: miss): reach missed target
% -23 = Undershoot (type 1: inaccurate): reach movement ended before target
% -24 = Undershoot (type 2: slow): reach was still ongoing when time expired, which occurred short of the target
% -31 = Overshoot (type 2: scuff): cursor barely scrapes the target and goes through it
% -32 = Overshoot (type 3: blow-through): reach goes right through the target
% -33 = Early Return: reaching back to center before time is up
% -34 = Target Hold Drift: hand drifts out of the end target
% -35 = Jitter: land at edge of target -> small jitter pulls cursor out
%
% -1x = delay epoch, -2x = reach epoch, -3x = target hold epoch
end

