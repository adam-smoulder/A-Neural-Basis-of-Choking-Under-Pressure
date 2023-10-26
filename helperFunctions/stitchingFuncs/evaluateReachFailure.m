function [trialStatusLabel] = evaluateReachFailure(curTime,goCueTime,failTime_wrtGo,...
    speed,endTargRadius,cursorRadius,rotPos,targetDistance,peakSpeed,peakSpeedInd,...
    postReachInd,reachEndInd,reachEndSpeed,reactionTimeInd,reactionTime,avgReachSpeedVal,varargin)
% For use in behaviorExtraction_forReach; evaluates what type of reach
% failure is occurring.
%
% Adam Smoulder, 3/10/20 (edit 6/18/20)

minPostReachSpeed = min(speed(peakSpeedInd:postReachInd));

% we use the minimum average speed as a way to see if
% reaches are demonstrably too slow.
%                 exitStartInd = find(distFromCenter > startTargRadius-cursorRadius, 1,'first'); % when the edge of the cursor leaves the center
%                 exitStartTime = curTime(exitStartInd);
%                 avgReachSpeed = mean(speed(reactionTimeInd:reachEndInd));
%                 minAvgReachSpeed = (distFromEndTarg(reactionTimeInd)-endTargRadius-startTargRadius+2*cursorRadius)/(failTime-exitStartTime);
minAvgReachSpeed = ((targetDistance-rotPos(reactionTimeInd,1))-endTargRadius)/(failTime_wrtGo-reactionTime); % assumes only target-axis distance needs covered and center of cursor in is good enough = conservative
% minAvgReachSpeed = (70-endTargRadius+cursorRadius)/(failTime-reactionTime); % RT is always achieved before 70mm distance unless it's a terrible reach, hence this is conservative

% To determine if it's a no attempt (that's well disguised in the form of a
% very lazy attempt), we need to know around when he exited the center
% target and if the peak speed occurred very late.
exitTime = curTime(find((sqrt(sum(rotPos.^2,2)) > 9) & curTime >= goCueTime,1));
peakSpeedTime = curTime(peakSpeedInd);

% To determine if the trial is a quit (reach in random direction), we need
% the angle to the location at time of failure (rotated so reach should be
% to 0 deg)
angleToFailRotPos = atan2d(rotPos(reachEndInd,2),rotPos(reachEndInd,1));

if peakSpeed < 0.2 || (((exitTime-goCueTime) < 400) && ((failTime_wrtGo+goCueTime-peakSpeedTime) < 100))
    % No Attempt: no attempt was made (or a very lazy one, indicated by a
    % super late peak speed even though the target was technically exited)
    trialStatusLabel = -21;
elseif abs(angleToFailRotPos) > 80
    % Quit-out: Intentional reach in wrong direction
    trialStatusLabel = -20;
elseif max(abs(rotPos(reactionTimeInd:reachEndInd,2))) > 25
    % Wild: Reach is very very not straight! 
    trialStatusLabel = -26;
elseif max(rotPos(peakSpeedInd:postReachInd,1)) > targetDistance
    % Inaccurate: reach missed the target
    trialStatusLabel = -22;
elseif minPostReachSpeed < 0.02 %  < 0.1*peakSpeed %
    % Inaccurate: main reach finished outside of the target
    trialStatusLabel = -23;
elseif avgReachSpeedVal < minAvgReachSpeed
    % Slow: Reach is too slow to have possibly made it
    trialStatusLabel = -24;
else
    % Slow: The remainder are on the ascent/descent of the
    % peak speed curve and are mid-reach.
    trialStatusLabel = -25;
end



% If other arguments are provided, plot it
if length(varargin)==6
    distFromEndTarg = varargin{1};
    startTargRadius = varargin{2};
    compareInds = varargin{3};
    ballisticPredVel = varargin{4};
    curBallistic5prcErrorTime = varargin{5};
    curRotatedBallisticPredEndpoint = varargin{6};
    
    
    reachEndTime = goCueTime+failTime_wrtGo;
    
    figure;
    subplot(3,2,[1 3])
    plot(curTime(1:postReachInd+18),distFromEndTarg(1:postReachInd+18),'.-','markersize',8)
    hold on
    % plot(curTime(exitStartInd), distFromEndTarg(exitStartInd),'b.','markersize',30)
    plot(curTime(reactionTimeInd), distFromEndTarg(reactionTimeInd),'r.','markersize',30)
    plot(curTime(peakSpeedInd), distFromEndTarg(peakSpeedInd),'m.','markersize',30)
    plot(curTime(reachEndInd), distFromEndTarg(reachEndInd),'k.','markersize',30)
    plot(curTime(compareInds(end)),sqrt(sum(([85 0]-curRotatedBallisticPredEndpoint).^2)),'k+','markersize',10,'color',[1,0.647,0],'linewidth',2)
    plot(goCueTime*ones(100,1),linspace(0,max(distFromEndTarg(1:postReachInd+18)),100),'g--')
    plot((failTime_wrtGo+goCueTime)*ones(100,1),linspace(0,max(distFromEndTarg(1:postReachInd+18)),100),'r--')
    plot(curTime(1:postReachInd+18),(endTargRadius)*ones(length(curTime(1:postReachInd+18)),1),'k--')
    xlabel('Time from Target Onset (ms)')
    axis([-inf inf -inf inf])
    ylabel('Distance from end target (mm)')
    title(['green = go, red = cutoff, black = target edge'])
    set(gca,'fontsize',16)
    legend('Data','Reaction Time','Peak Speed','Reach End (or end of trial)','Ball. Pred. Endpt','Go Cue','Fail Time','End Targ Radius','location','SW')
    
    
    
    subplot(3,2,[2 4])
    plot(curTime(1:postReachInd+18),speed(1:postReachInd+18))
    hold on
    plot(reachEndTime*ones(100,1),linspace(min(speed(1:postReachInd+18)),max(speed(1:postReachInd+18)),100),'r--')
    plot(goCueTime*ones(100,1),linspace(min(speed(1:postReachInd+18)),max(speed(1:postReachInd+18)),100),'g--')
    % plot(curTime(exitStartInd), speed(exitStartInd),'b.','markersize',30)
    plot(curTime(reactionTimeInd), speed(reactionTimeInd),'r.','markersize',30)
    plot(curTime(peakSpeedInd), peakSpeed,'m.','markersize',30)
    plot(curTime(reachEndInd), reachEndSpeed,'k.','markersize',30)
    plot(curTime(compareInds),sqrt(sum(ballisticPredVel.^2,2)),'k--','linewidth',1,'color',[1,0.647,0])
    plot(curBallistic5prcErrorTime,0,'ko','markersize',10,'color',[1,0.647,0],'linewidth',2)
    axis([-inf inf -inf inf])
    ylabel('Speed (m/s)')
    % title(['s_{avg,min} = ' num2str(round(minAvgReachSpeed,2)), ...
    %     ', s_{avg} = ' num2str(round(avgReachSpeed,2)),...
    %     ', s_{peak} = ' num2str(round(peakSpeed,2)),...
    %     ', s_{minPostReach} = ' num2str(round(minPostReachSpeed,2))])
    infoString1 = ['RT = ' num2str(reactionTime) 'ms,  ',...
        's_{avg,min} = ' num2str(round(minAvgReachSpeed,2)), ...
        ',  s_{avg} = ' num2str(round(avgReachSpeedVal,2))];
    infoString2 = ['s_{peak} = ' num2str(round(peakSpeed,2)),...
        ',  s_{minPostReach} = ' num2str(round(minPostReachSpeed,2))];
    text(50+curTime(1), 0.9*peakSpeed,infoString1)
    text(50+curTime(1), 0.8*peakSpeed,infoString2)
    
    set(gca,'fontsize',16)
    
    
    
    subplot(3,2,[5 6])
    plot(rotPos((1:postReachInd+18),1),rotPos((1:postReachInd+18),2),'k.-','markeredgecolor','b','markersize',10)
    theta = linspace(0,2*pi,100);
    startTarget = startTargRadius*[cos(theta); sin(theta)]';
    target = endTargRadius*[cos(theta); sin(theta)]'+[85 0];
    cursorFailInd = find(curTime < reachEndTime,1,'last')+1;
    cursor = cursorRadius*[cos(theta); sin(theta)]'+rotPos(cursorFailInd,:);
    hold on
    plot(startTarget(:,1),startTarget(:,2),'g-')
    plot(target(:,1),target(:,2),'k-')
    plot(cursor(:,1),cursor(:,2),'r-')
    % plot(rotPos(exitStartInd,1),rotPos(exitStartInd,2),'b.','markersize',30)
    plot(rotPos(reactionTimeInd,1),rotPos(reactionTimeInd,2),'r.','markersize',30)
    plot(rotPos(peakSpeedInd,1),rotPos(peakSpeedInd,2),'m.','markersize',30)
    plot(rotPos(reachEndInd,1),rotPos(reachEndInd,2),'k.','markersize',30)
    plot(curRotatedBallisticPredEndpoint(1),curRotatedBallisticPredEndpoint(2),'k+','markersize',10,'color',[1,0.647,0],'linewidth',2)
    set(gca,'fontsize',16)
    xlabel('Target Axis (mm)')
    ylabel('Off-Target Axis (mm)')
    title(['Class = ' num2str(trialStatusLabel)])
    
%     set(gcf,'position',[1921 41 1920 963]) % for lab setup
    set(gcf,'position',[1 41 1680 933]) % for WFH setup
end

end

