function [lEstO,masterInfo,behavior] = realignDataPostStitching(lEstO,masterInfo,behavior,alignOption,roundToTimeBin)
% The purpose of this function is to realign the data in a specific way.
% This typically assumes the data uses small time bins (e.g. 5ms), though
% it can technically be done with larger ones. 
%
% For all alignment options, it first see if the trial in masterInfo exists
% in behavior. If it does, it sees if it has the desired alignment point
% (e.g. go cue). If it has this point (it is either a successful reach or a
% trial that at least got to this point), that will be used for alignment.
% If not, it will use the time of failure to align.
%
% Inputs: 
%
% - lEstO: [ntimebins x nfactors] matrix of orthonormalized factor
% scores. We will remove values for stuff that can't be aligned.
%
% - masterInfo: [ntimebins x nlabels] Matrix with many different labels.
% A few matter here specifically: The 1st column must be the day label,
% the 2nd column must be the trial label, the 3rd column should be the
% state label, and the 4th column should be the time label (within
% trial). Because we'll be aligning based on 'behavior', we do not assume
% any given temporal alignment of the data (see below). WE ASSUME MASTER
% INFO IS DOUBLES, not int16 (sometimes used for storage)
%
% - behavior: [ntrial x 1] structure with all relevant behavioral stuff.
% There's a lot. Not gonna list it all here. We're mainly going to be
% taking advantage of the field on behavior that start with "t#", "day",
% "trial", and "trialStatusLabels"
%
% - alignOption: One of the different alignment options. There are a lot,
% and all are listed below. Should be a string.
%
% - roundToTimeBin: boolean; if true, alignment will be performed such that
% some value in masterInfo(:,4) is 0 for each trial, where 0 is the closest
% bin to the desired alignment point. If false, the alignments will just be
% subtracted for each trial, meaning that there will not be a 0 in the time
% column of masterInfo unless the exact value of a neural timebin is the
% same as the desired alignment.
%
% 
% Outputs:
% - lEstO: same as input but with unalignable stuff removed
% - masterInfo: same as input but each trial has been realigned
% - behavior: same as input but each trial has been realigned
%
%
% Choices for alignState:
%
% - 'delay' or 'targetOnset': aligns to when the target shows up. This
% should be present on all trials (t2_targetOnsetTimes)
%
% - 'goCueOrFail': for trials that make it to go cue, it aligns to that.
% For trials that don't, it aligns to the time of failure. This is
% t3_goCueTimes.
%
% - 'reachStartOrFail': for trials that make it to reach, aligns based on
% t4b_reactionTimes_40mm. NOTE THAT THIS IS A MISNOMER it's actually
% 20mm/s as the speed. Otherwise, align to time of failure
%
% - 'peakSpeedOrFail': time of peak speed for trials that have a reach.
% Otherwise time of failure.
%
% - 'targetEntryOrFail': time of target entry for trials that have a reach.
% Otherwise time of failure. Uses t7_postReachStateTimes
%
% - 'reachEndOrFail': if there is a full reach, uses
% t6b_reachEndTimes_40mm (similarly misnomeric as explained above).
% Otherwise uses failure time.
%
% - 'successTimeOrFail': if success, aligns to end of hold period
% (t8_trialEndTimes). Otherwise, to time of failure.
%
%
% Adam Smoulder, 4/14/20 (edit 4/15/20)


% If "delay" is the choice, rename it to "targetOnset"
if strcmp(alignOption,'delay'), alignOption = 'targetOnset'; end

% We first get stuff off masterInfo and behavior that we'll need
dayTrialLabels_factor = masterInfo(:,1)*10000+masterInfo(:,2);
dayTrials = unique(dayTrialLabels_factor);
timeLabels_factor = masterInfo(:,4);
dayTrialLabels_behavior = double([behavior.day])*10000+double([behavior.trial]);

% We loop over each dayTrial in the data and align based on behavior. If
% there is no alignment, we add it to indsToRemove and get rid of it.
indsToRemove = [];
for i = 1:length(dayTrials)
    curDayTrial = dayTrials(i);
    curInds = find(dayTrialLabels_factor==curDayTrial);
    curBehInd = find(dayTrialLabels_behavior==curDayTrial,1,'first');
    if isempty(curBehInd) % missing; can't align off of it
        indsToRemove = [indsToRemove curInds];
    else % we got it; let's try to align
        switch alignOption % pick field to align to
            case 'targetOnset'
                timeToCheckName = 't2_targetOnsetTimes';
            case 'goCueOrFail'
                timeToCheckName = 't3_goCueTimes';
            case 'reachStartOrFail'
                timeToCheckName = 't4b_reactionTimes_40mm';
            case 'peakSpeedOrFail'
                timeToCheckName = 't5_peakSpeedTimes';
            case 'reachEndOrFail'
                timeToCheckName = 't6b_reachEndTimes_40mm';
            case 'targetEntryOrFail'
                timeToCheckName = 't7_postReachStateTimes';
            case 'successTimeOrFail'
                timeToCheckName = 't8_trialEndTimes';
            otherwise
                disp('Unusable align option selected')
                timeToCheckName = [];
        end
        
        alignTime = [];
        if isempty(timeToCheckName) % if we choose something that doesn't exist (in this trial), skip
            indsToRemove = [indsToRemove curInds];
        else % otherwise, align it!
            % if the value is nan, it's a failure; we take the failure time
            if isnan(behavior(curBehInd).(timeToCheckName))
                if behavior(curBehInd).trialStatusLabels > -30 % if fail before end of reach
                    timeToCheckName = 't7_postReachStateTimes'; % When failure is triggered for fail trials
                else % failure after reach stopped
                    timeToCheckName = 't8_trialEndTimes';
                end
            end
            
            % Next, if we need to round to the nearest masterInfo time bin, we do so.
            if roundToTimeBin
                curTime = timeLabels_factor(curInds); % get masterInfo times for the current indices
                [~,closestValInd] = min(abs(curTime-double(behavior(curBehInd).(timeToCheckName)))); % find index of closest value
                alignTime = curTime(closestValInd); % set the align time to this value
            else
                alignTime = double(behavior(curBehInd).(timeToCheckName));
            end
            
            % apply to masterInfo
            masterInfo(curInds,4) = masterInfo(curInds,4)-alignTime;
            
            % since behavior is in integers mostly, we cast alignTime to
            % that and use whichever is needed by field
            alignTimeInt = int16(round(alignTime));
            alignTimeDouble = double(alignTime);
            
            % apply to all the behavior fields 
            behavior(curBehInd).t1_trialStartTimes = behavior(curBehInd).t1_trialStartTimes-alignTimeDouble;
            behavior(curBehInd).t2_targetOnsetTimes = behavior(curBehInd).t2_targetOnsetTimes-alignTimeDouble;
            behavior(curBehInd).t3_goCueTimes = behavior(curBehInd).t3_goCueTimes-alignTimeDouble;
            behavior(curBehInd).t4_reactionTimes = behavior(curBehInd).t4_reactionTimes-alignTimeDouble;
            behavior(curBehInd).t5_peakSpeedTimes = behavior(curBehInd).t5_peakSpeedTimes-alignTimeDouble;
            behavior(curBehInd).t6_reachEndTimes = behavior(curBehInd).t6_reachEndTimes-alignTimeDouble;
            behavior(curBehInd).t7_postReachStateTimes = behavior(curBehInd).t7_postReachStateTimes-alignTimeDouble;
            behavior(curBehInd).t8_trialEndTimes = behavior(curBehInd).t8_trialEndTimes-alignTimeDouble;
            behavior(curBehInd).kinematics_updated.time =  behavior(curBehInd).kinematics_updated.time-alignTimeDouble;
            tempStateTrans = behavior(curBehInd).stateTrans;
            if iscell(tempStateTrans)
                behavior(curBehInd).stateTrans{1}(2,:) = tempStateTrans{1}(2,:)-alignTimeDouble;
            else
                behavior(curBehInd).stateTrans(2,:) = tempStateTrans(2,:)-alignTimeDouble;
            end
            
            if isfield(behavior(curBehInd),'t4b_reactionTimes_40mm')
                behavior(curBehInd).t4b_reactionTimes_40mm = behavior(curBehInd).t4b_reactionTimes_40mm-alignTimeDouble;
            end
            if isfield(behavior(curBehInd),'t4b_reactionTimes_20mm')
                behavior(curBehInd).t4b_reactionTimes_20mm = behavior(curBehInd).t4b_reactionTimes_20mm-alignTimeDouble;
            end
            if isfield(behavior(curBehInd),'t4c_reactionTimes_exit')
                behavior(curBehInd).t4c_reactionTimes_exit = behavior(curBehInd).t4c_reactionTimes_exit-alignTimeDouble;
            end
            if isfield(behavior(curBehInd),'t6b_reachEndTimes_40mm')
                behavior(curBehInd).t6b_reachEndTimes_40mm = behavior(curBehInd).t6b_reachEndTimes_40mm-alignTimeDouble;
            end
            if isfield(behavior(curBehInd),'t6b_reachEndTimes_20mm')
                behavior(curBehInd).t6b_reachEndTimes_20mm = behavior(curBehInd).t6b_reachEndTimes_20mm-alignTimeDouble;
            end
            if isfield(behavior(curBehInd),'eyeData_updated') && isfield(behavior(curBehInd),'eyeDataIsGood')
                if behavior(curBehInd).eyeDataIsGood==1
                    behavior(curBehInd).eyeData_updated.time = behavior(curBehInd).eyeData_updated.time-alignTimeDouble;
                end
            end
        end
    end
    
    % display progress sometimes
    if mod(i,500)==0
        disp(['Realigned trial ' num2str(i) ' of ' num2str(length(dayTrials))])
    end
end; clear i

% Remove invalid indices
lEstO(indsToRemove,:) = [];
masterInfo(indsToRemove,:) = [];

disp(['Completed realignment to ' alignOption ' ' ternop(roundToTimeBin,'with','without') ' rounding to neural time bins'])
end

