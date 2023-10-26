function [lines] = plotAcrossTimeBins(data,timebins)
% plotAcrossTimeBins makes plots that show a given metric across both delay
% and reach period epochs. From the input, it infers that there are a total
% of T timebins and R desired traces. 
%
% This DOES both create a new figure and use hold on; change this if you
% feel...
%
% Timebin setup is as follows:
% The first 0 represents target onset, while the 2nd 0 represents the start 
% of the first reach bin.  For instance:
% [-200 0 200 400 -400 -200 0 200 400] would show the bin before
% target onset thru 2 bins after target onset, then 2 bins before reach
% onset thru 2 bins after reach onset.
%
% READ HERE: If the data has a 3rd dimension (N), these are assumed to be
% repetitions from a BOOTSTRAPPED distribution (cuz that's the data I'm
% normally looking at here...). Shading will represent the SEM
% of this.
%
% Inputs:
% - data:  [T x R x N] data, where N is assumed to be bootstrapped samples
% of each of the R traces along timebins T. The mean across N will be
% plotted for each of the R traces and will be shaded by the SEM over N
% (empirical) at each timebin.
% - timebins:  [T x 1] list of timebins for both the delay and reach
% period. The transition between the two is assumed to be when the 1st
% difference goes negative for a bin.
%
% Outputs:
% - lines: variables for mainTraces such that you can set legends to them

% get info, traces, CI, other info for plotting
[T,R,N] = size(data);
t = 1:T; % we'll relabel later
mainTraces = squeeze(mean(data,3));
semTraces = squeeze(std(data,[],3)); % It's SEM bc the data is bootstrapped
%[lowTraces, highTraces] = bootConfInt(data,95);

% find important timepoints and axis values
targetOnsetInd = find(timebins==0,1,'first');
reachOnsetInd = find(timebins==0,1,'last');
transitionInd = find(diff(timebins) < 0, 1,'first'); % this is the last target onset centered bin
maxVal = max(mainTraces(:)+semTraces(:));
minVal = min(mainTraces(:)-semTraces(:));

colors = getDistinctColors;
lines = [];

figure
hold on
plot(targetOnsetInd*ones(1,100),linspace(minVal,maxVal,100),'r--','linewidth',3)
plot(reachOnsetInd*ones(1,100),linspace(minVal,maxVal,100),'g--','linewidth',3)
for i = 1:R
    curLine = plot(mainTraces(:,i),'-','color',colors{i},'linewidth',3);
    lines = [lines ; curLine];
    if N > 1
        %CITraces = ([highTraces(:,i)' ; lowTraces(:,i)']-mainTraces(:,i)').*[1; -1];
        shadedErrorBar(t,mainTraces(:,i),semTraces(:,i),'lineProps',{'r-','color',colors{i}},'linewidth',1.5)
    end
end
plot((transitionInd+0.5)*ones(1,100),linspace(minVal,maxVal,100),'k-','linewidth',10)
xticks(t)
xticklabels(num2str(timebins))
axis([-inf inf minVal maxVal])
set(gca,'fontsize',18)
if N > 1
    title('SEM error shading')
end
xlabel('Time from event (ms)')

end

