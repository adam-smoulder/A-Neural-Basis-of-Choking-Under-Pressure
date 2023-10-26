function [trialStatusLabel] = evaluateDelayFailure(distFromCenter,startTargRadius,rotPos,postReachStateTime,trialStartTime,curTime,targetOnsetTime,varargin)
% For use in behaviorExtraction_forReach; evaluates what type of delay
% failure is occurring.
%
% Adam Smoulder, 3/10/20 (edit 3/12/20)

failInd = find(curTime > postReachStateTime,1,'first');
periFailInds = (failInd-6):(failInd+18);

if (max(distFromCenter(periFailInds)) > 1.75*startTargRadius) && (max(rotPos(periFailInds,1)) > startTargRadius) && (postReachStateTime-targetOnsetTime) > 100
    % False start: attempted reach towards target early
    trialStatusLabel = -11;
elseif (max(distFromCenter(periFailInds)) <= 1.75*startTargRadius) && (postReachStateTime-targetOnsetTime) > 100
    % Cheat/drift: no full reach attempt made in any dir
    trialStatusLabel = -12;
elseif (postReachStateTime-targetOnsetTime) <= 100
    % Mis-start: accidentally enters/exits start target
    trialStatusLabel = -13;
elseif max(rotPos(periFailInds,1)) < 1/4*startTargRadius
    % Quit out: reach made away from target
    trialStatusLabel = -14;
else
    % Shouldn't happen; if it does, need to readjust criteria
    trialStatusLabel = -111;
end



% If other information is given, we can plot it
if length(varargin)==1
    speed = varargin{1};
%     targetOnsetTime = varargin{2};
    lastIndToPlot = periFailInds(end);
    trialStartInd = find(curTime > trialStartTime,1,'first');
    targetOnsetInd = find(curTime > targetOnsetTime,1,'first');
    
    figure;
    subplot(3,2,[1 3])
    plot(curTime(trialStartInd:lastIndToPlot),distFromCenter(trialStartInd:lastIndToPlot),'.-','linewidth',2,'markersize',8)
    hold on
    plot(targetOnsetTime*ones(100,1),linspace(min(distFromCenter(trialStartInd:lastIndToPlot)),max(distFromCenter(trialStartInd:lastIndToPlot)),100),'m--')
    plot(postReachStateTime*ones(100,1),linspace(min(distFromCenter(trialStartInd:lastIndToPlot)),max(distFromCenter(trialStartInd:lastIndToPlot)),100),'r--')
    plot(curTime(trialStartInd:lastIndToPlot),startTargRadius*ones(length((trialStartInd:lastIndToPlot)),1),'k--')
    axis([-inf inf -inf inf])
    title(['Distance from start targ (mm)'])
    set(gca,'fontsize',16)
    
    
    subplot(3,2,[2 4])
    plot(curTime(trialStartInd:lastIndToPlot),speed(trialStartInd:lastIndToPlot),'linewidth',2)
    hold on
    plot(targetOnsetTime*ones(100,1),linspace(min(speed(trialStartInd:lastIndToPlot)),max(speed(trialStartInd:lastIndToPlot)),100),'m--')
    plot(postReachStateTime*ones(100,1),linspace(min(speed(trialStartInd:lastIndToPlot)),max(speed(trialStartInd:lastIndToPlot)),100),'r--')
    axis([-inf inf -inf inf])
    title('Speed (m/s)')
    set(gca,'fontsize',16)
    
    
    subplot(3,2,[5 6])
    plot(rotPos((trialStartInd:lastIndToPlot),1),rotPos((trialStartInd:lastIndToPlot),2),'k.-','markeredgecolor','b','markersize',10)
    theta = linspace(0,2*pi,100);
    startTarget = startTargRadius*[cos(theta); sin(theta)]';
    target = 7*[cos(theta); sin(theta)]'+[85 0];
    cursorFailInd = find(curTime < postReachStateTime,1,'last')+1;
    cursor = 3*[cos(theta); sin(theta)]'+rotPos(cursorFailInd,:);
    hold on
    plot(startTarget(:,1),startTarget(:,2),'g-')
    plot(target(:,1),target(:,2),'k-')
    plot(cursor(:,1),cursor(:,2),'r-')
    plot(rotPos(targetOnsetInd,1),rotPos(targetOnsetInd,2),'m.','markersize',30)
    plot(rotPos(cursorFailInd,1),rotPos(cursorFailInd,2),'r.','markersize',30)
    set(gca,'fontsize',16)
    xlabel('Target Axis (mm)')
    ylabel('Off-Target Axis (mm)')
    title(['Class = ' num2str(trialStatusLabel)])
    
%     set(gcf,'position',[1921 41 1920 963]) % for lab setup
    set(gcf,'position',[1 41 1680 933]) % for WFH setup
    
end

end

