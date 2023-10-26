function [trialData] = updateReactionTimeB(trialData,taskInfo)
% This is a post-processing hack to improve our estimate of reaction time;
% originally one of the reaction times was the time when the hand speed was
% 0.02 m/s; however, this isn't super robust to quick movements around.
% Making it 0.02 m/s towards the target seems a good bit fairer, which is
% what we do here.
%
% Inputs: 
% - trialData: [ntrials x 1] structure that comes out of stitching
% preprocessing.
% - taskInfo: structure from stitching preprocessing; has the angles
% associated with the direction labels
%
% Output:
% - trialData: [ntrials x 1] same as input, but now the field 
% "t4b_reactionTime_20mm" is updated to be the time it takes to achieve a
% hand speed of 20mm/s in the direction of the target instead of period
%
% Adam Smoulder, 2/9/21

threshSpeed = 0.02; % 0.02 m/s = 20 mm/s

for i = 1:length(trialData)
    curTrial = trialData(i);
    goCueTime = curTrial.t3_goCueTime;
    if ~isnan(goCueTime)
        % Get stuff we need, like peak speed time and speed towards the
        % target
        time = curTrial.time;
        peakSpeedTime = curTrial.t5_peakSpeedTime;
        dir = curTrial.directionLabel;
        curVel = curTrial.kinematics_updated.velocity;
        curAngle = taskInfo.directionLabelAngles(dir);
        curRotMat = [cosd(curAngle) -sind(curAngle) ; sind(curAngle) cosd(curAngle)];
        curRotVel = curVel*curRotMat;
        speedTowardsTarget = curRotVel(:,1);
        
        % Going backward from peak speed, find when the speedTowardsTarget
        % first exceeds 0.02 m/s; this time is the new RTb
        indsToCheck = time <= peakSpeedTime & time >= goCueTime;
        timeToCheck = time(indsToCheck);
        RTInd = find(speedTowardsTarget(indsToCheck)<threshSpeed,1,'last')+1;
        if isempty(RTInd)
            trialData(i).t4b_reactionTime_20mm = nan;
        else
            trialData(i).t4b_reactionTime_20mm = timeToCheck(RTInd);
        end
    end
end; clear i

end

