function [trialStatusLabel] = evaluateTargetHoldFailure(curTime,reachEndTime,...
    stateTrans,endTimeInd,rotPos,speed,endTargRadius,...
    targetHoldEndInd,enterInd,timeInTarget,distInTarget,targDist,varargin)
% For use in behaviorExtraction_forReach; evaluates what type of reach
% failure is occurring.
%
% Adam Smoulder, 3/10/20 (edit 3/12/20)

% hold failure
theta = linspace(0,2*pi,100);
target = endTargRadius*[cos(theta); sin(theta)]'+[targDist 0];
exitSpeed = mean(speed(targetHoldEndInd:min((targetHoldEndInd+5),length(curTime))));
% exitVec = diff(rotPos(failInd-1:failInd,:)); % + = away from center, - = towards center
% exitDir = atan2(exitVec(2),exitVec(1));

[~,closestCircleInd] = min(sum(((target-rotPos(targetHoldEndInd,:)).^2),2)); % minimize distance from exit to target pt
exitDir = theta(closestCircleInd);
if exitDir > pi % we'll use the negative convention
    exitDir = -(2*pi-exitDir);
end

if isempty(enterInd)
    % This is a bug...it never really entered, should happen rarely. Call
    % it a slow reach.
    trialStatusLabel = -333;
elseif targetHoldEndInd-enterInd <= 1
    % This tends to be a system bug...for some reason, the trial is a
    % failure immediately on entry, even if it doesn't really look like a
    % scuff. Not sure why. Note it and throw it out.
    trialStatusLabel = -334;
elseif (timeInTarget <= 125) && (distInTarget < endTargRadius)
    % Scuff: Cursor just barely scrapes the target
    trialStatusLabel = -31;
elseif (timeInTarget < 250) && (distInTarget > endTargRadius) && (abs(exitDir) < 3*pi/4)
    % Overshoot: Reach goes straight through target
    trialStatusLabel = -32;
elseif (exitSpeed > 0.1) && (abs(exitDir) > pi/2)
    % Early Return: Reaching back to center before time is up
    trialStatusLabel = -33;
elseif distInTarget > endTargRadius
    % Drift: Doesn't blow through target, but doesn't hold properly, drifting out of the target
    trialStatusLabel = -34;
else
    % Jitter: Holding at edge of target -> jitters out
    trialStatusLabel = -35;
end




% plot if we have additional argument (true)
if length(varargin) > 0
    compareInds = varargin{1};
    ballisticPredVel = varargin{2};
    curBallistic5prcErrorTime = varargin{3};
    curRotatedBallisticPredEndpoint = varargin{4};
    
    % (curTime,reachEndTime,stateTrans,endTimeInd,rotPos,speed,endTargRadius)
    enterInd = find(curTime > reachEndTime,1,'first');
    distFromEndTarg = sqrt(sum((rotPos-[targDist 0]).^2,2));
    lastIndToPlot = find(curTime<(reachEndTime+400),1,'last')+1;
    if lastIndToPlot > length(curTime)
        lastIndToPlot = targetHoldEndInd;
    end
    
    figure;
    subplot(3,2,[1 3])
    plot(curTime(1:lastIndToPlot),distFromEndTarg(1:lastIndToPlot),'.-','markersize',8)
    hold on
    plot(curTime(targetHoldEndInd), distFromEndTarg(targetHoldEndInd),'k.','markersize',30)
    plot(reachEndTime*ones(1000,1),linspace(0,max(distFromEndTarg),1000),'r--')
    plot((reachEndTime+400)*ones(1000,1),linspace(0,max(distFromEndTarg),1000),'b--')
    plot(linspace(curTime(1),curTime(lastIndToPlot),1000),(endTargRadius)*ones(1000,1),'k--')
    plot(curTime(compareInds(end)),sqrt(sum(([targDist 0]-curRotatedBallisticPredEndpoint).^2)),'k+','markersize',10,'color',[1,0.647,0],'linewidth',2)
    xlabel('Time from Target Onset (ms)')
    axis([-inf inf -inf inf])
    ylabel('Distance from end target (mm)')
    title(['Red = entry, Black = fail, blue = when success would have been'])
    set(gca,'fontsize',16)
    legend('Data','End of Trial','Target Entry','Reward Time','R.Targ Radius','Pred. Ball. Endpt')
    % legend('Data','Reaction Time','Peak Speed','Reach End (or end of trial)','location','SW')

    subplot(3,2,[2 4])
    plot(curTime(1:lastIndToPlot),speed(1:lastIndToPlot))
    hold on
    plot(reachEndTime*ones(100,1),linspace(min(speed(1:targetHoldEndInd)),max(speed(1:targetHoldEndInd)),100),'r--')
    plot((reachEndTime+400)*ones(100,1),linspace(min(speed(1:targetHoldEndInd)),max(speed(1:targetHoldEndInd)),100),'b--')
    plot(curTime(targetHoldEndInd), speed(targetHoldEndInd),'k.','markersize',30)
    plot(curTime(compareInds),sqrt(sum(ballisticPredVel.^2,2)),'k--','linewidth',1,'color',[1,0.647,0])
    plot(curBallistic5prcErrorTime,0,'ko','markersize',10,'color',[1,0.647,0],'linewidth',2)
    axis([-inf inf -inf inf])
    ylabel('Speed (m/s)')
    infoString1 = ['Time in target = ' num2str(timeInTarget) 'ms,  ',...
        'distance in target = ' num2str(distInTarget)];
    infoString2 = ['Exit Speed = ' num2str(exitSpeed),...
        ', exit direction = ' num2str(exitDir)];
%     text(50, 0.9*max(speed(1:targetHoldEndInd)),infoString1)
%     text(50, 0.8*max(speed(1:targetHoldEndInd)),infoString2)
    set(gca,'fontsize',16)
    
    subplot(3,2,[5 6])
    plot(rotPos((1:lastIndToPlot),1),rotPos((1:lastIndToPlot),2),'k.-','markeredgecolor','b','markersize',10)
    % startTarget = startTargRadius*[cos(theta); sin(theta)]';
    startTarget = 8.5*[cos(theta); sin(theta)]';
    cursor = 3*[cos(theta); sin(theta)]'+rotPos(targetHoldEndInd,:);
    hold on
    plot(startTarget(:,1),startTarget(:,2),'g-')
    plot(target(:,1),target(:,2),'k-')
    plot(cursor(:,1),cursor(:,2),'r-')
    plot(rotPos(enterInd,1),rotPos(enterInd,2),'r.','markersize',30)
    plot(rotPos(targetHoldEndInd,1),rotPos(targetHoldEndInd,2),'k.','markersize',30)
    plot(curRotatedBallisticPredEndpoint(1),curRotatedBallisticPredEndpoint(2),'k+','markersize',10,'color',[1,0.647,0],'linewidth',2)
    set(gca,'fontsize',16)
    xlabel('Target Axis (mm)')
    ylabel('Off-Target Axis (mm)')
    title(['Class = ' num2str(trialStatusLabel)])
    
    %     set(gcf,'position',[1921 41 1920 963]) % for lab setup
    set(gcf,'position',[1 41 1680 933]) % for WFH setup
    disp('plotting shit')
end

