function [behavior] = behaviorExtraction(data,trialStartState,goCueState,reachState,...
    targHoldState,failTime,center,startTargRadius,endTargRadius,...
    possibleRewardSizes,targetDistance,targetOnPostReactionState,...,
    cursorRadius,allDistPrcs,dayNum)
% This function runs behavior extraction. I don't like using it as a
% function instead of a script, but it's the only damn way to keep track of
% it in all the different contexts and stuff I use it.
%
% This is meant to be used after standard preprocessing.
%
% Adam Smoulder, 6/9/20

theta = linspace(0,2*pi,100);
rotatedTarget = endTargRadius*[cos(theta); sin(theta)]'+[targetDistance 0];

% Things to save
stateTrans = cell(1);          % state transitions; captures everything timing-wise
distPrctileTimes = cell(1); % vector of times when the animal reaches [5:5:95]% of the distance to the target...
% We save this in a 3 x nprcs matrix, where the top row is allDistPrcs, the middle row is the distances, and the bottom row is the times

% Get alll the important timepoints for the reach. the "t#" preceding means
% it's w.r.t. the trial and not a raw quantity (e.g. reactionTime here is
% not true reaction time; it would be RT = t4_reactionTimes-t3_goCueTimes)
t1_trialStartTimes = [];       % what is the first timepoint for this trial?
t2_targetOnsetTimes = [];      % when did the target turn on?
t3_goCueTimes = [];            % when did the go cue happen?
t4_reactionTimes = [];         % when was 20% of peak speed achieved?
t4b_reactionTimes_20mm = [];   % when did speed rise to 0.02m/s = 20 mm/s?
t4c_reactionTimes_exit = [];   % when did the cursor exit the center target?
t5_peakSpeedTimes = [];        % when was peak speed achieved?
t6_reachEndTimes = [];         % when was 10% of peak speed achieved (after peak speed)? OR when timeout occurs
t6b_reachEndTimes_20mm = [];   % when did speed fall to 0.02m/s = 20 mm/s (roughly 5% of avg peak speed)
t7_postReachStateTimes = [];   % when did the state change from reach to either failure or target hold?
t8_trialEndTimes = [];         % when did reward/failure occur?

bSpeed = 0.02; % speed to use for reaction time and reach end time B

% other important single metrics
centerExitSpeed = [];       % what was the velocity upon exiting the center target?
kinematics_updated = [];    % kinematics after correcting for center and inverted y-axis + other things
eyeDataIsGood = [];         % Boolean flag on if the eye data for the given trial is usable or not
eyeData_updated = [];       % eye data after correcting for center and inverted y-axis + other things
peakSpeed = [];             % what was the maximum speed in the reach?
avgReachSpeed = [];         % what was the overall average speed during the reach?
reactionTime = [];          % how quickly did the animal reach 20% of peak speed after the go cue?
timeInTarget = [];          % how long was the cursor in the end target?
distInTarget = [];          % how much distance did the cursor cover in the target?
trialStatusLabels = [];     % what was the outcome of the trial? All possible outcomes are listed above



% Alright let's actually get all of this stuff...
counter = 0;
rewLengths = []; % we don't save this, but use it sometimes for quick checks
dirs = [];
days = [];
trials = [];

% Get targets and reorder them
allTargs = nan(min(length(data),500),2); % get all possible reach targets
for j = 1:length(allTargs)
    j
    allTargs(j,:) = data(j).TrialData.endTarget(1:2);
end; clear j
reachTargs = unique(allTargs,'rows');
[~,sortOrder] = sort(reachTargs(:,1),'descend');
reachTargs = reachTargs(sortOrder,:); % highest X should be target 1
temp = reachTargs-center;
targAngles = round(atan2d(-temp(:,2),temp(:,1)));
targAngles(targAngles < 0) = targAngles(targAngles < 0) + 360;
if size(reachTargs,1)==8
    anglesShouldBe = (0:45:315)';
elseif size(reachTargs,1)==2
    anglesShouldBe = [0;180];
end
[~,order] = ismember(anglesShouldBe,targAngles);
reachTargs = reachTargs(order,:)
    
% Now actually get stuff for each trial
for j = 1:length(data)
    if contains(data(j).Overview.trialName,'Choice Task') % we only want single target
        continue
    end
    counter = counter+1;
    
    % Get state transitions; determine if unattempted; also get direction
    % and reward labels
    curStateTrans = data(j).TrialData.stateTransitions;
    stateTrans{counter} = curStateTrans;
    if isfield(data(j),'codes') % assignments preproc has been done
        rewLengths(counter) = data(j).codes(2);
        dirs(counter) = data(j).codes(1);
    else % for behavior paper, get reward and direction labels
        states = data(j).Parameters.stateNames;
        stateMatrix = nan(length(states),1);
        stateMatrix = num2cell(stateMatrix);
        for k = 1:size(stateMatrix,1)
            stateMatrix{k,1} = states{1,k};
        end
        targetHoldIndex = find(strcmp('Target Hold',stateMatrix) == 1);
        if isempty(targetHoldIndex)
            targetHoldIndex = targHoldState;
        end
        rewardPassState = data(j).Parameters.StateTable(targetHoldIndex).windowPassState;
        if iscell(possibleRewardSizes(1)) % using names, found in data(j).Definitions.AllPossibleIntervals.names
            curRewName = data(j).Parameters.StateTable(rewardPassState).Interval.name{1};
            rewLengths(counter) = find(strcmp(possibleRewardSizes,curRewName)); % 1 is the first reward, 2 is the second, etc.
        else % is a number = using length of reward
            rewardSize = data(j).Parameters.StateTable(rewardPassState).Interval.length;   % Reward size of Reach
            rewLengths(counter) = find(possibleRewardSizes==rewardSize);
        end
        
        
        %             [~,~,dirlabel] = intersect(data(j).TrialData.endTarget(1:2),reachTargs,'rows');
        %             dirs(counter) = dirlabel;
        % HACK for when the center changes..URIOSUEFGHUIOSDG
        temp = data(j).TrialData.endTarget(1:2)-center;
        targAngle = round(atan2d(-temp(:,2),temp(:,1)));
        targAngle(targAngle < 0) = targAngle(targAngle < 0) + 360;
        [~,dirlabel] = min(abs(anglesShouldBe-targAngle));
        dirs(counter) = dirlabel;
    end
    days(counter) = dayNum;
    trials(counter) = str2double(data(j).Overview.trialNumber(6:end));
    
    
    if ~any(curStateTrans(1,:)==targetOnPostReactionState) % reach never shown; call unattempted
        trialStatusLabels(counter) = 0;
        t1_trialStartTimes(counter) = nan;
        t2_targetOnsetTimes(counter) = nan;
        t3_goCueTimes(counter) = nan;
        t4_reactionTimes(counter) = nan;
        t4b_reactionTimes_20mm(counter) = nan;
        t4c_reactionTimes_exit(counter) = nan;
        t5_peakSpeedTimes(counter) = nan;
        t6_reachEndTimes(counter) = nan;
        t6b_reachEndTimes_20mm(counter) = nan;
        t7_postReachStateTimes(counter) = nan;
        t8_trialEndTimes(counter) = nan;
        distPrctileTimes(counter) = {nan};
        centerExitSpeed(counter) = nan;
        peakSpeed(counter) = nan;
        avgReachSpeed(counter) = nan;
        reactionTime(counter) = nan;
        timeInTarget(counter) = nan;
        distInTarget(counter) = nan;
        kinematics_updated(counter).position = nan;
        kinematics_updated(counter).velocity = nan;
        kinematics_updated(counter).acceleration = nan;
        kinematics_updated(counter).time = nan;
        kinematics_updated(counter).rotatedPosition = nan;
        kinematics_updated(counter).distanceFromCenter = nan;
        kinematics_updated(counter).distanceFromEndTarget = nan;
        kinematics_updated(counter).speed = nan;
        eyeDataIsGood(counter) = false;
        eyeData_updated(counter).position = nan;
        eyeData_updated(counter).velocity = nan;
        eyeData_updated(counter).time = nan;
        eyeData_updated(counter).rotatedPosition = nan;
        eyeData_updated(counter).distanceFromCenter = nan;
        eyeData_updated(counter).distanceFromEndTarget = nan;
        eyeData_updated(counter).speed = nan;
        disp(['Unattempted trial, day ' num2str(j) ', trial ' num2str(j)])
        continue
    end
    
    % If attempted, let's find when different trial checkpoints happen
    trialStartStateInd = find(curStateTrans(1,:)==trialStartState,1,'first');
    targetOnsetStateInd = find(curStateTrans(1,:)==targetOnPostReactionState-1,1,'first');
    goCueStateInd = find(curStateTrans(1,:)==goCueState,1,'first');
    postReachStateInd = find(curStateTrans(1,:)>reachState,1,'first');
    if isempty(trialStartStateInd) % bug
        t1_trialStartTimes(counter) = nan;
    else % normal
        t1_trialStartTimes(counter) = double(curStateTrans(2,trialStartStateInd));
    end % unattempted; probably shouldn't get to here...
    if isempty(targetOnsetStateInd)
        t2_targetOnsetTimes(counter) = nan;
    else
        t2_targetOnsetTimes(counter) = double(curStateTrans(2,targetOnsetStateInd));
    end
    if isempty(goCueStateInd) % fail before go cue
        t3_goCueTimes(counter) = nan;
    else % normal
        t3_goCueTimes(counter) = double(curStateTrans(2,goCueStateInd));
    end
    if isempty(postReachStateInd) % this should always exist?
        t7_postReachStateTimes(counter) = nan;
    else
        t7_postReachStateTimes(counter) = double(curStateTrans(2,postReachStateInd));
    end
    
    % Get kinematic stuff we need
    j
    kins = data(j).TrialData.HandKinematics.markerKinematics; % NOTE for some data its capital H on HandKinematics...ugh
    indsToUse = kins.time > (curStateTrans(2,trialStartState)-9); % about one sample before trial start state
    curTime = kins.time(indsToUse);
    curPos = (kins.position(indsToUse,1:2)-center).*[1 -1];
    curVel = kins.velocity(indsToUse,1:2).*[1,-1];
    curAcc = kins.acceleration(indsToUse,1:2).*[1,-1];
    targLoc = (data(j).TrialData.endTarget(1:2)-center).*[1 -1]; % y-axes are flipped for mirror
    targVec = targLoc/sum(sqrt(targLoc.^2));
    targAngle = atan2d(targVec(2),targVec(1));
    rotMat = [cosd(targAngle) -sind(targAngle) ; sind(targAngle) cosd(targAngle)];
    curRotPos = curPos*rotMat;
    curDistFromCenter = sqrt(curPos(:,1).^2+curPos(:,2).^2);
    curDistFromEndTarg = sqrt((targLoc(1)-curPos(:,1)).^2+(targLoc(2)-curPos(:,2)).^2);
    curSpeed = sqrt(kins.velocity(indsToUse,1).^2+kins.velocity(indsToUse,2).^2);
    
    % Get eye data stuff we want
    if isfield(data(j).TrialData,'EyeData')
        eyeData = data(j).TrialData.EyeData;
        eyeData.position = eyeData.position;
        % We're going to chop off the first and last time index to
        % better estimate velocity for the same timebins
        if ~isempty(eyeData.position) && length(eyeData.position(:,1))~=1
            eyeDataIsGood(counter) = true;
            eyeData_updated(counter).position = (eyeData.position(2:end-1,1:2)-center).*[1 -1];
            rawVel = (diff(eyeData.position(:,1:2)).*[1 -1])./diff(eyeData.time);
            eyeData_updated(counter).velocity = (rawVel(1:end-1,:)+rawVel(2:end,:))/2;
            eyeData_updated(counter).time = eyeData.time(2:end-1);
            eyeData_updated(counter).rotatedPosition = eyeData_updated(counter).position*rotMat;
            eyeData_updated(counter).distanceFromCenter = sqrt(eyeData_updated(counter).position(:,1).^2+eyeData_updated(counter).position(:,2).^2);
            eyeData_updated(counter).distanceFromEndTarget = sqrt((targLoc(1)-eyeData_updated(counter).position(:,1)).^2+(targLoc(2)-eyeData_updated(counter).position(:,2)).^2);
            eyeData_updated(counter).speed = sqrt(eyeData_updated(counter).velocity(:,1).^2+eyeData_updated(counter).velocity(:,2).^2);
        else
            eyeDataIsGood(counter) = false;
        end
        if ~isempty(eyeData.pupilSize) && sum(isnan(eyeData.pupilSize))==0
            eyeData_updated(counter).pupilSize = eyeData.pupilSize(2:end-1);
        end
    else
        eyeDataIsGood(counter) = false;
        eyeData_updated(counter).position = nan;
        eyeData_updated(counter).velocity = nan;
        eyeData_updated(counter).time = nan;
        eyeData_updated(counter).rotatedPosition = nan;
        eyeData_updated(counter).distanceFromCenter = nan;
        eyeData_updated(counter).distanceFromEndTarget = nan;
        eyeData_updated(counter).speed = nan;
        eyeData_updated(counter).pupilSize = nan;
    end
    % Get indices for events w.r.t. curTime
    %       trialStartInd = find(curTime >= t1_trialStartTimes(end),1,'first');
    targetOnsetInd = find(curTime >= t2_targetOnsetTimes(end),1,'first');
    %       goCueInd = find(curTime >= t3_goCueTimes(end),1,'first');
    postReachInd = find(curTime >= t7_postReachStateTimes(end),1,'first');
    
    % Get exit speed stuff
    upperRTExitInd = find(curDistFromCenter > startTargRadius & curTime >= t2_targetOnsetTimes(end),1,'first');
    lowerRTExitInd = upperRTExitInd-1;
    if isempty(upperRTExitInd) || (lowerRTExitInd==0)
        t4c_reactionTimes_exit(counter) = nan;
        centerExitSpeed(counter) = nan;
    else
        t4c_reactionTimes_exit(counter) = twoPointInterpolation(startTargRadius,...
            curDistFromCenter(lowerRTExitInd),curDistFromCenter(upperRTExitInd),...
            curTime(lowerRTExitInd),curTime(upperRTExitInd));
        centerExitSpeed(counter) = twoPointInterpolation(startTargRadius,...
            curDistFromCenter(lowerRTExitInd),curDistFromCenter(upperRTExitInd),...
            curSpeed(lowerRTExitInd),curSpeed(upperRTExitInd));
    end
    % Get reach speed / timing stuff
    [curPeakSpeed, peakSpeedInd] = max(curSpeed(targetOnsetInd:postReachInd));
    peakSpeedInd = peakSpeedInd+targetOnsetInd-1; % readjust to be inline with curTime and everything
    peakSpeed(counter) = curPeakSpeed;
    reactionTimeInd = find(curSpeed(1:peakSpeedInd) < 0.2*curPeakSpeed,1,'last')+1;
    if isempty(reactionTimeInd) % probably means it's crap...
        reactionTimeInd = find(curSpeed(1:peakSpeedInd) >= 0.2*curPeakSpeed,1,'first');
    end
    reactionTime(counter) = curTime(reactionTimeInd)-t3_goCueTimes(end);
    lowerRT20Ind = find(curSpeed(1:peakSpeedInd) < bSpeed,1,'last');
    upperRT20Ind = lowerRT20Ind+1;
    if isempty(lowerRT20Ind)
        t4b_reactionTimes_20mm(counter) = nan;
    else
        t4b_reactionTimes_20mm(counter) = twoPointInterpolation(bSpeed,curSpeed(lowerRT20Ind),curSpeed(upperRT20Ind),curTime(lowerRT20Ind),curTime(upperRT20Ind));
    end
    %         reactionTime40mmInd = find(curSpeed(1:peakSpeedInd) < 0.04,1,'last')
    %         t4b_reactionTimes40mm(counter) = curTime(reactionTime40mmInd);
    reachStateDuration = t7_postReachStateTimes(end)-t3_goCueTimes(end); % important for determining failure mode
    t4_reactionTimes(counter) = curTime(reactionTimeInd);
    t5_peakSpeedTimes(counter) = curTime(peakSpeedInd);
    if ismember(targHoldState,curStateTrans(1,:)) % if end target hold was achieved
        reachEndInd = find((curSpeed < 0.1*curPeakSpeed) ...
            & (1:length(curSpeed) > peakSpeedInd)'...
            ,1,'first'); % slowing down to 10% of peak speed = reach end
        if isempty(reachEndInd) % odd... rarely occurs. Just use the last time I supposed
            reachEndInd = length(curSpeed);
        end
    else % a reach or targ hold failure; cap it at the maximum reach time
        reachEndInd = find((curSpeed < 0.1*curPeakSpeed) ...
            & (1:length(curSpeed) > peakSpeedInd)'...
            & (1:length(curSpeed) <= postReachInd)'...
            ,1,'first');
        if isempty(reachEndInd)
            reachEndInd = postReachInd;
        end
    end
    
    t6_reachEndTimes(counter) = curTime(reachEndInd);
    avgReachSpeed(counter) = mean(curSpeed(reactionTimeInd:reachEndInd));
    
    upperReachEnd20Ind = find((curSpeed < bSpeed) ...
        & (1:length(curSpeed) > peakSpeedInd)'...
        ,1,'first'); % slowing down to 10% of peak speed = reach end
    lowerReachEnd20Ind = upperReachEnd20Ind-1;
    if isempty(upperReachEnd20Ind)
        t6b_reachEndTimes_20mm(counter) = nan;
    else
        t6b_reachEndTimes_20mm(counter) = twoPointInterpolation(bSpeed,curSpeed(lowerReachEnd20Ind),curSpeed(upperReachEnd20Ind),curTime(lowerReachEnd20Ind),curTime(upperReachEnd20Ind));
    end
    
    reachEndSpeed = curSpeed(reachEndInd);
    
    
    % Get how long the cursor was in the target and how much it moved
    % during the end target hold interval
    if ismember(targHoldState,curStateTrans(1,:)) % if end target hold was achieved
        % Find when success/failure occurs (end of target hold)
        targetHoldEndInd = find(curTime > curStateTrans(2,postReachStateInd+1), 1,'first');
        if isempty(targetHoldEndInd) % on some successes, the reward is the last timepoint; hence we don't have any kinematics (curTime) following
            targetHoldEndInd = length(curTime);
        end
        
        % Find the entry point by searching the 5 pts previous (and 4 pts later) to timer start
        valsToCheck = curRotPos(postReachInd-4:postReachInd+4,:); % get relevant vals
        
        % calculate if the target-axis value is above that needed to be in target for the given height
        enterInd = find(valsToCheck(:,1)>(85-sqrt(endTargRadius.^2-valsToCheck(:,2).^2)),1,'first')-5+postReachInd;
        if isempty(enterInd) % if it is empty, this is an error; we correct for it later
            timeInTarget(counter) = nan;
            distInTarget(counter) = nan;
        else
            timeInTarget(counter) = curTime(targetHoldEndInd)-curTime(enterInd-1);
            distInTarget(counter) = sum(sqrt(sum(diff(curRotPos(enterInd:targetHoldEndInd,:)).^2,2)));
        end
    else
        timeInTarget(counter) = nan;
        distInTarget(counter) = nan;
    end
    
    % Calculate the times when different percentiles of the distance
    % gap are covered. We look at percentile of ideal distance to
    % target.
    distToCover = targetDistance; % don't bother correcting; do so later. %-(startTargRadius+endTargRadius)+cursorRadius;
    curDistPrcTimes = nan(length(allDistPrcs),1);
    allDistForPrcs = nan(length(allDistPrcs),1);
    for nprc = 1:length(allDistPrcs)
        tempPrcDist = (1-allDistPrcs(nprc)/100)*distToCover; % + endTargRadius; % the distance we're looking for
        allDistForPrcs(nprc) = tempPrcDist;
        upperDistPrcInd = find(curTime >= t2_targetOnsetTimes(end) & ...
            curDistFromEndTarg <= tempPrcDist,1,'first');
        lowerDistPrcInd = upperDistPrcInd-1;
        if ~isempty(upperDistPrcInd) && lowerDistPrcInd ~= 0
            curDistPrcTimes(nprc) = twoPointInterpolation(tempPrcDist,...
                curDistFromEndTarg(upperDistPrcInd),curDistFromEndTarg(lowerDistPrcInd),...
                curTime(upperDistPrcInd),curTime(lowerDistPrcInd));
        end
    end; clear nprc
    distPrctileTimes(counter) = {[allDistPrcs ; allDistForPrcs' ; curDistPrcTimes']};
    
    
    % We calculate a predicted ballistic reach endpoint. We do this by
    % mirroring the velocity traces up to peak speed and finding the
    % endpoint. We also want to calculate the "time of deviation" from
    % this mirrored trajectory; we'll use error as 10% of peak speed.
    
    % First, we mirror the trajectory; we do this from when the speed
    % exceeds 0.02ms to peak speed
    % offset = 6; % roughly 3 indices = 25ms
    %         if isempty(lowerRT20Ind)
    %             mirrorStartInd = reactionTimeInd; % only should happen for delay failures...
    %         else
    %             mirrorStartInd = lowerRT20Ind;
    %         end
    if isempty(peakSpeedInd) || (peakSpeedInd-30)<1
        mirrorStartInd = reactionTimeInd; % only should happen for delay failures...
    else
        mirrorStartInd = peakSpeedInd-30;
    end
    mirrorInds = mirrorStartInd:peakSpeedInd; % mirrorInds = (reactionTimeInd-offset):peakSpeedInd;
    velToMirror = curVel(mirrorInds,:);
    ballisticPredVel = [velToMirror(1:end-1,:) ; flip(velToMirror)]; % don't double count peak speed
    compareInds = mirrorStartInd:mirrorStartInd-1+length(ballisticPredVel); % reactionTimeInd-offset:reactionTimeInd-offset-1+length(ballisticPredVel);
    if max(compareInds)>length(curVel)
        compareInds = compareInds(1:find(compareInds==length(curVel),1,'first'));
        ballisticPredVel = ballisticPredVel(1:find(compareInds==length(curVel),1,'first'),:);
    end
    realVelToCompare = curVel(compareInds,:); % the actual velocity at this time
    ballisticError_normPS = abs(realVelToCompare-ballisticPredVel)/peakSpeed(end);
    
    % Find where error exceeds 5%
    upper5percInd = find(sqrt(sum(ballisticError_normPS.^2,2)) > 0.05,1,'first');
    lower5percInd = upper5percInd-1;
    if isempty(upper5percInd)
        curBallistic5prcErrorTime = curTime(compareInds(end)+1);
    elseif lower5percInd == 0 % bad trial...
        curBallistic5prcErrorTime = 0;
    else
        curBallistic5prcErrorTime = twoPointInterpolation(0.05,...
            ballisticError_normPS(lower5percInd),ballisticError_normPS(upper5percInd),...
            curTime(compareInds(lower5percInd)),curTime(compareInds(upper5percInd)));
    end
    
    % Find where error exceeds 10%
    upper10percInd = find(sqrt(sum(ballisticError_normPS.^2,2)) > 0.1,1,'first');
    lower10percInd = upper5percInd-1;
    if isempty(upper10percInd)
        curBallistic10prcErrorTime = curTime(compareInds(end)+1);
    elseif lower10percInd == 0 % bad trial...
        curBallistic10prcErrorTime = 0;
    else
        curBallistic10prcErrorTime = twoPointInterpolation(0.1,...
            ballisticError_normPS(lower10percInd),ballisticError_normPS(upper10percInd),...
            curTime(compareInds(lower10percInd)),curTime(compareInds(upper10percInd)));
    end
    
    % Find where error exceeds 20%
    upper20percInd = find(sqrt(sum(ballisticError_normPS.^2,2)) > 0.2,1,'first');
    lower20percInd = upper5percInd-1;
    if isempty(upper20percInd)
        curBallistic20prcErrorTime = curTime(compareInds(end)+1);
    elseif lower20percInd == 0 % bad trial...
        curBallistic20prcErrorTime = 0;
    else
        curBallistic20prcErrorTime = twoPointInterpolation(0.2,...
            ballisticError_normPS(lower20percInd),ballisticError_normPS(upper20percInd),...
            curTime(compareInds(lower20percInd)),curTime(compareInds(upper20percInd)));
    end
    
    % Calculate the predicted endpoint by mirroring displacement similarly
    dispToPS = curPos(peakSpeedInd,:)-curPos(mirrorInds(1),:);
    curBallisticPredEndpoint = dispToPS*2+curPos(mirrorInds(1),:);
    curRotatedBallisticPredEndpoint = curBallisticPredEndpoint*rotMat;
    
    % assign the trial status based on the type of failure/success
    if double(data(j).Overview.trialStatus)==1 % success
        trialStatusLabels(counter) = 1;
        t8_trialEndTimes(counter) = double(curStateTrans(2,postReachStateInd+1));
%                           evaluateTargetHoldFailure(curTime,t7_postReachStateTimes(end),curStateTrans,postReachStateInd,curRotPos,curSpeed,endTargRadius,targetHoldEndInd,enterInd,timeInTarget(end),distInTarget(end),targetDistance,compareInds,ballisticPredVel,curBallistic20prcErrorTime,curRotatedBallisticPredEndpoint); % to plot; same setup as target hold fail
%                           subplot(3,2,[5 6]); title(['Class = 1 dirs = ' num2str(dirs(end)) ' rews = ' num2str(rewLengths(end))]) % but some of the titles are messed up
%                           disp('yep')
        
    else % is a failure
        % delay failure
        if isnan(reachStateDuration) || ((reachStateDuration < failTime-100) && ~any(ismember(curStateTrans(1,:),targHoldState))) % -100 is such that any possible little muckups in timing don't falsely say "oh yeah, reach duration was < failTime" when it's just because of 1/2s resolution errors
            trialStatusLabels(counter) = evaluateDelayFailure(curDistFromCenter,startTargRadius,curRotPos,t7_postReachStateTimes(end),t1_trialStartTimes(end),curTime,t2_targetOnsetTimes(end));
            %                                 evaluateDelayFailure(curDistFromCenter,startTargRadius,curRotPos,t7_postReachStateTimes(end),t1_trialStartTimes(end),curTime,t2_targetOnsetTimes(end),curSpeed); % used for plotting
            %                                 disp(['Day ' num2str(i) ' trial ' num2str(j)])
            %                                 disp(['dirs = ' num2str(dirs(end)) ' rews = ' num2str(rewLengths(end))])
            %                                 disp('yep')
            t3_goCueTimes(end) = nan;
            t4_reactionTimes(end) = nan;
            t4b_reactionTimes_20mm(end) = nan;
            t5_peakSpeedTimes(end) = nan;
            t6_reachEndTimes(end) = nan;
            t6b_reachEndTimes_20mm(end) = nan;
            peakSpeed(end) = nan;
            avgReachSpeed(end) = nan;
            reactionTime(end) = nan;
            t8_trialEndTimes(counter) = double(curStateTrans(2,postReachStateInd));
            if ismember(trialStatusLabels(end),[-111]) % plot weirdos...
%                 evaluateDelayFailure(curDistFromCenter,startTargRadius,curRotPos,t7_postReachStateTimes(end),t1_trialStartTimes(end),curTime,curSpeed,t2_targetOnsetTimes(end)); % used for plotting
                disp(['-111 Bad Trial - day ' num2str(j) ' trial ' num2str(j)])
            end
            
            % reach fail
        elseif (abs(reachStateDuration-failTime)<2) && ~any(ismember(curStateTrans(1,:),targHoldState)) % it's a time out in the reach
            trialStatusLabels(counter) = evaluateReachFailure(curTime,t3_goCueTimes(end),failTime,curSpeed,endTargRadius,cursorRadius,curRotPos,targetDistance,curPeakSpeed,peakSpeedInd,postReachInd,reachEndInd,reachEndSpeed,reactionTimeInd,reactionTime(end),avgReachSpeed(end));
            t8_trialEndTimes(counter) = double(curStateTrans(2,postReachStateInd));
            %                                 evaluateReachFailure(curTime,t3_goCueTimes(end),failTime,curSpeed,endTargRadius,cursorRadius,curRotPos,targetDistance,curPeakSpeed,peakSpeedInd,postReachInd,reachEndInd,reachEndSpeed,reactionTimeInd,reactionTime(end),avgReachSpeed(end),curDistFromEndTarg,startTargRadius,compareInds,ballisticPredVel,curBallistic20prcErrorTime,curRotatedBallisticPredEndpoint); % for plotting
            %                                 disp(['Day ' num2str(i) ' trial ' num2str(j)])
            %                                 disp(['dirs = ' num2str(dirs(end)) ' rews = ' num2str(rewLengths(end))])
            %                                 disp('yep')
            
            
            % hold failure
        elseif any(ismember(curStateTrans(1,:),targHoldState))
            trialStatusLabels(counter) = evaluateTargetHoldFailure(curTime,t7_postReachStateTimes(end),curStateTrans,postReachStateInd,curRotPos,curSpeed,endTargRadius,targetHoldEndInd,enterInd,timeInTarget(end),distInTarget(end),targetDistance);
            t8_trialEndTimes(counter) = double(curStateTrans(2,postReachStateInd+1));
            if trialStatusLabels(end)==-333 % magic number = this is actually a reach failure
                disp(['Falsely "landed" reach; these are jitters, day ' num2str(j) ' trial ' num2str(j)])
                %                     trialStatusLabels(end) = evaluateReachFailure(curTime,goCueTimes(end),failTime,speed,endTargRadius,cursorRadius,rotPos,targetDistance);
                trialStatusLabels(end) = -35;
                t8_trialEndTimes(counter) = double(curStateTrans(2,postReachStateInd));
            elseif trialStatusLabels(end)==-334 % Magic number - throw this one away
                % These trials tend to look like scuffs or jitters, but
                % it's such a short period of time in the target that
                % it's hard to tell.
                disp(['Bad trial (odd failure?), throwing away, day ' num2str(j) ' trial ' num2str(j)])
                trialStatusLabels(end) = 0;
            end
%             if rewLengths(end)==4
%                                             evaluateTargetHoldFailure(curTime,t7_postReachStateTimes(end),curStateTrans,postReachStateInd,curRotPos,curSpeed,endTargRadius,targetHoldEndInd,enterInd,timeInTarget(end),distInTarget(end),targetDistance,compareInds,ballisticPredVel,curBallistic20prcErrorTime,curRotatedBallisticPredEndpoint);
%                                             disp(['Day ' num2str(i) ' trial ' num2str(j)])
%                                             disp(['dirs = ' num2str(dirs(end)) ' rews = ' num2str(rewLengths(end))])
%                                             disp('yep')
%             end
            
        else
            trialStatusLabels(counter) = 0; % ??? This shouldn't happen, we should've caught it earlier
        end
    end
    
    
    if curBallistic5prcErrorTime-t5_peakSpeedTimes(end) > 260
        disp('Sample drops for ballistic error calc led to >260ms; calling 260')
        disp(['Day ' num2str(j) ' trial ' num2str(j)])
        curBallistic5prcErrorTime = t5_peakSpeedTimes(end)+260
    end
    
    if curBallistic10prcErrorTime-t5_peakSpeedTimes(end) > 260
        disp('Sample drops for ballistic error calc led to >260ms; calling 260')
        disp(['Day ' num2str(j) ' trial ' num2str(j)])
        curBallistic10prcErrorTime = t5_peakSpeedTimes(end)+260
    end
    
    if curBallistic20prcErrorTime-t5_peakSpeedTimes(end) > 260
        disp('Sample drops for ballistic error calc led to >260ms; calling 260')
        disp(['Day ' num2str(j) ' trial ' num2str(j)])
        curBallistic20prcErrorTime = t5_peakSpeedTimes(end)+260
    end
    
    % make new kinematics structure for ease of use
    newKins.position = curPos;
    newKins.velocity = curVel;
    newKins.acceleration = curAcc;
    newKins.time = curTime;
    
    % Update the state transitions and kinematics stuff
    kinematics_updated(counter).position = curPos;
    kinematics_updated(counter).velocity = curVel;
    kinematics_updated(counter).acceleration = curAcc;
    kinematics_updated(counter).time = curTime;
    kinematics_updated(counter).rotatedPosition = curRotPos;
    kinematics_updated(counter).distanceFromCenter = curDistFromCenter;
    kinematics_updated(counter).distanceFromEndTarget = curDistFromEndTarg;
    kinematics_updated(counter).speed = curSpeed;
    kinematics_updated(counter).ballistic5prcErrorTime = curBallistic5prcErrorTime-t5_peakSpeedTimes(end); % correcting
    kinematics_updated(counter).ballistic10prcErrorTime = curBallistic5prcErrorTime-t5_peakSpeedTimes(end); % correcting
    kinematics_updated(counter).ballistic20prcErrorTime = curBallistic20prcErrorTime-t5_peakSpeedTimes(end); % correcting
    kinematics_updated(counter).ballisticPredEndpoint = curBallisticPredEndpoint;
    kinematics_updated(counter).rotatedBallisticPredEndpoint = curRotatedBallisticPredEndpoint;
end; clear j



% Assign fields to behavior
task_parameters.center = center;
task_parameters.cursorRadius = cursorRadius;
task_parameters.startTargRadius = startTargRadius;
task_parameters.endTargRadius = endTargRadius;
task_parameters.targetDistance = targetDistance;
task_parameters.trialStartState = trialStartState;
task_parameters.targetOnPostReactionState = targetOnPostReactionState;
task_parameters.goCueState = goCueState;
task_parameters.reachState = reachState;
task_parameters.targHoldState = targHoldState;
task_parameters.failTime = failTime;

for j = 1:length(trialStatusLabels)
    behavior(j).direction = dirs(j);
    behavior(j).reward = rewLengths(j);
    
    behavior(j).trialStatusLabels = trialStatusLabels(j);      % what was the outcome of the trial? All possible outcomes are:
    behavior(j).reactionTimes = reactionTime(j);               % how quickly did the animal reach 20% of peak speed after the go cue?
    behavior(j).peakSpeed = peakSpeed(j);                      % what was the maximum speed in the reach?
    behavior(j).avgReachSpeed = avgReachSpeed(j);              % what was the overall average speed during the reach?
    behavior(j).centerExitSpeed = centerExitSpeed(j);
    
    behavior(j).t1_trialStartTimes = t1_trialStartTimes(j);            % what is the first timepoint for this trial?
    behavior(j).t2_targetOnsetTimes = t2_targetOnsetTimes(j);          % when did the target turn on?
    behavior(j).t3_goCueTimes = t3_goCueTimes(j);                      % when did the go cue happen?
    behavior(j).t4_reactionTimes = t4_reactionTimes(j);                % when was 20% of peak speed achieved?
    behavior(j).t4b_reactionTimes_20mm = t4b_reactionTimes_20mm(j);    % when did speed rise to 0.04m/s = 40 mm/s?
    behavior(j).t4c_reactionTimes_exit = t4c_reactionTimes_exit(j);
    behavior(j).t5_peakSpeedTimes = t5_peakSpeedTimes(j);              % when was peak speed achieved?
    behavior(j).t6_reachEndTimes = t6_reachEndTimes(j);                % when was 10% of peak speed achieved (after peak speed)? OR when timeout occurs
    behavior(j).t6b_reachEndTimes_20mm = t6b_reachEndTimes_20mm(j);    % when did speed fall to 0.04m/s = 40 mm/s (roughly 10% of avg peak speed)
    behavior(j).t7_postReachStateTimes = t7_postReachStateTimes(j);    % when did the state change from reach to either failure or target hold?
    behavior(j).t8_trialEndTimes = t8_trialEndTimes(j);                % when did reward/failure occur?
    
    behavior(j).kinematics_updated = kinematics_updated(j);     % updated kinematics after recentering and correcting y-axis inversion
    behavior(j).eyeDataIsGood = eyeDataIsGood(j);               % is the eye data for this trial usable? boolean
    behavior(j).eyeData_updated = eyeData_updated(j);           % updated eye data after recentering and correcting y-axis inversion
    behavior(j).stateTrans = stateTrans{j};                     % state transitions = ; captures everything timing-wise
    behavior(j).distPrctileTimes = distPrctileTimes{j};
    behavior(j).distInTarget = distInTarget(j);                 % how much distance did the cursor cover in the target?
    behavior(j).timeInTarget = timeInTarget(j);                 % how long was the cursor in the end target?
    
    behavior(j).task_parameters = task_parameters;
    behavior(j).trialIndex = trials(j);
end; clear j


end

