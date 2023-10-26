function behavior = behaviorExtraction_preStitching(fnames,animalName)
% This function is used in assignConditions_Earl.m and extracts relevant
% behavioral info from the data, assigning it into the output "behavior"
% variable. This includes various timepoints in the task, various kinematic
% metrics (e.g. reaction time, peak speed), an updated kinematics
% structure, and detailed trial status labels (e.g. overshoot, freeze/quit)
%
% Inputs:
% - fnames: cell array of preprocessed data file names
% - monkeyName: used to determine which parameter set to use (e.g. target
% size)
%
% Outputs:
% - behavior: contains relevant behavioral metrics and info. See bottom of
% function for all fields
%
% Adam Smoulder, 3/16/20 (edit 3/24/20)
%
% Detailed trial status labels can take on the following values:
%%% 1 = success
%%% 0 = unattempted
%%% -11 = False start: attempted reach towards target early
%%% -12 = Cheat/drift: no full reach attempt made in any dir
%%% -13 = Mis-start: accidentally enters/exits start target
%%% -14 = Quit out: reach made away from target
%%% -21 = Freeze/quit: no attempt was made
%%% -22 = Inaccurate: reach missed the target (overshot and missed)
%%% -23 = Inaccurate: main reach finished outside of the target (misplace but no overshoot)
%%% -24 = Slow: Reach is too slow to have possibly made it
%%% -25 = Slow: The remainder are on the ascent/descent of the peak speed curve and are mid-reach.
%%% -31 = Scuff: Cursor just barely scrapes the target
%%% -32 = Overshoot: Reach goes straight through target
%%% -33 = Early Return: Reaching back to center before time is up
%%% -34 = Drift: Doesn't blow through target, but doesn't hold properly, drifting out of the target
%%% -35 = Jitter: Holding at edge of target -> jitters out


if strcmp(animalName,'Nelson')
    % Nelson
    trialStartState = 2;
    goCueState = 5;
    reachState = 6;
    targHoldState = 7;
    failTime = 625;
%     center = [120 370];
%     startTargRadius = 8.7;
%     endTargRadius = 5.7;
    possibleRewardSizes = [75 122 170 980];
    targetDistance = 85;
    
elseif strcmp(animalName,'Ford')
    
    % Ford, no tubes
    trialStartState = 2;
    goCueState = 5;
    reachState = 6;
    targHoldState = 7;
    failTime = 625;
%     center = [-10 400];
%     startTargRadius = 8.35;
%     endTargRadius = 5.25;
    possibleRewardSizes = [84 135 177 901];
    targetDistance = 85;
    
elseif strcmp(animalName,'EarlB')
    % Earl
    trialStartState = 2;
    goCueState = 5;
    reachState = 6;
    targHoldState = 7;
    failTime = 750;
%     center = [90,455];
%     startTargRadius = 8.3;
%     endTargRadius = 7.3;
    possibleRewardSizes = [84 135 177 901]; % only used for OLD (for behavior paper)
    targetDistance = 85;
    
elseif strcmp(animalName,'EarlN') || strcmp(animalName,'Earl')
    trialStartState = 2;
    goCueState = 5;
    reachState = 6;
    targHoldState = 7;
    failTime = 750;
%     center = [90,455];
%     startTargRadius = 8.3;
%     endTargRadius = 7.3;
    possibleRewardSizes = [0 135 221 901]; % use for NEW earl data (neural paper)
    targetDistance = 85; 
    
elseif strcmp(animalName,'EarlB_centerOut') % Use 20200608-20200609
    trialStartState = 2;
    goCueState = 3;
    reachState = 3;
    targHoldState = 4;
    failTime = 5000;
%     center = [90,430];
%     startTargRadius = 8;
%     endTargRadius = 7;
    possibleRewardSizes = [160];
    targetDistance = 90;
end

% consistent across all tasks/animals
targetOnPostReactionState = 4;
cursorRadius = 3;
allDistPrcs = 10:1:100; % percentiles of distance to target to calc time at

for i = 1:length(fnames)
    load([pwd '\' fnames{i}])
    startTargRadius = data(1).TrialData.startTarget(end);
    center = data(1).TrialData.startTarget(1:2);
    endTargRadius = data(1).TrialData.endTarget(end);
    dayNum = i;
    dayBehavior = behaviorExtraction(data,trialStartState,goCueState,reachState,...
        targHoldState,failTime,center,startTargRadius,endTargRadius,...
        possibleRewardSizes,targetDistance,targetOnPostReactionState,...
        cursorRadius,allDistPrcs,dayNum);
    if i == 1
        behavior = dayBehavior;
    else
        behavior = [behavior dayBehavior];
    end
    disp(['Go cue times and trial statuses acquired for day ' num2str(i)])
end; clear i data
        
end

