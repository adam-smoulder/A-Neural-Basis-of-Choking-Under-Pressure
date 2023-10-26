function [kinematicsAligned] = alignKinematicsToMasterInfo(kinematics,masterInfo)
% This function interpolates 2D kinematics from the timestamps given to
% timestamps aligning to the masterInfo. 
%
% Inputs:
% - kinematics: ndays x 1 cell array, where each cell contains an array of
% ntrials x 1 structures. The structure for each given trial has 4
% properties: time, position, velocity, and acceleration. Each of these
% fields is ntimebins x 1.
% - masterInfo: (total) nbins x 9 matrix, where nbins represents every
% single time bin recorded across all sessions. The columns of master info
% are as follows: 
% 1) day #, 2) trial #, 3) state #, 4) time in trial, 5) tube #, 6) reward
% #, 7) punishment #, 8) Success/fail = 1/0, 9) single vs. 2 target task =
% 0/1.
% In this function, we only use day, trial, and time.
%
% Outputs:
% - kinematicsAligned: structure with same 4 properties as before, though
% now each is nbins x 1 and matches with the time column of masterInfo
%
% I'm sure there's a much more efficient way to write all of this,
% but it takes a (survivable) ~3s per day. So run when you have a
% half hour to kill and save the result.
%
% Adam Smoulder, 9/17/19

nbins = size(masterInfo,2);
ndays = length(kinematics);
dayLabels = masterInfo(:,1);
trialLabels = masterInfo(:,2);
timeLabels = masterInfo(:,4);

kinTimeAligned = nan(nbins,1);
kinPosAligned =  nan(nbins,3);
kinVelAligned =  nan(nbins,3);
kinAccAligned =  nan(nbins,3);

% frustratingly, I didn't save any info on the trial # with the kinematics;
% this means we'll have to check that there's the same number of trials in
% the masterInfo as there is in kinematics and assume they're all aligned
% (they should be). If there's a diff number...er...we'll have to go back
%tic
bincounter = 1; % to keep track of where we are w.r.t. masterInfo
for d = 1:ndays
    curDayKin = kinematics{d};
    ntrials = length(curDayKin);
    trials = unique(trialLabels(dayLabels==d));
    assert(length(trials)==ntrials); % this is the assertion I mention above
    tempTimeAligned = []; % we do this weird temp method because it makes stuff wayyyyy faster. 
    tempPosAligned = [];  % I guess there are some sorta memory pressure problems with constantly assigning values
    tempVelAligned = [];  % to huge matrices (e.g. nbins x 3), as using these "temporary" variables that are cleared
    tempAccAligned = [];  % each day for assignment (so size nbins/ndays x 3) takes it from 90s/day down to 3-4s/day.
    oldbincounter = bincounter;
    for j = 1:ntrials
        % get kinematic stuff for the trial
        curKin = curDayKin(j); % sometimes this line has to be changed from curly to normal parentheses or vice versa
        if iscell(curKin), curKin = curKin{1}; end
        curKinTime = curKin.time;
        curKinPos = curKin.position;
        curKinVel = curKin.velocity;
        curKinAcc = curKin.acceleration;
        
        % get masterInfo stuff for the trial
        curInds = (dayLabels==d) & (trialLabels==trials(j));
        curTimeLabels = timeLabels(curInds);
        
        % there will always be more kinematics than masterInfo unless the
        % bin size for masterInfo is < screen Hz. Hence, we'll find the two
        % kinematic points around each masterInfo time bin and
        % weighted-average them. Same as linear interpolation of 2 pts.
        for t = 1:length(curTimeLabels)
            curTime = curTimeLabels(t);
            mapInd = find(curKinTime >= curTime,1,'first');
            if isempty(mapInd) % use the last one then; we're at the end
                mapInd = length(curKinTime);
            end
            if mapInd==1 % use the first bin with a lil hack
                mapInd = 2;
                weights = [1 0]; % so we just use the first bin; 0 weight for the 2nd!
            else
                weights = [abs(curKinTime(mapInd-1)-curTime) abs(curKinTime(mapInd)-curTime)];
            end
            weights = weights./sum(weights);
            tempTimeAligned(end+1) = curTime;
            tempPosAligned(end+1,:) = [curKinPos(mapInd-1,:) ; curKinPos(mapInd,:)]' * weights'; % weighted avg; [3 2] * [1 2]'
            tempVelAligned(end+1,:) = [curKinVel(mapInd-1,:) ; curKinVel(mapInd,:)]' * weights';
            tempAccAligned(end+1,:) = [curKinAcc(mapInd-1,:) ; curKinAcc(mapInd,:)]' * weights';
           
            bincounter = bincounter+1;
        end; clear t
    end; clear j
    kinTimeAligned(oldbincounter:bincounter-1) = tempTimeAligned;
    kinPosAligned(oldbincounter:bincounter-1,:) = tempPosAligned;
    kinVelAligned(oldbincounter:bincounter-1,:) = tempVelAligned;
    kinAccAligned(oldbincounter:bincounter-1,:) = tempAccAligned;
    clear tempTimeAligned tempVelAligned tempAccAligned tempPosAligned
    disp(['Aligned kinematics for day ' num2str(d)])
    %toc
end; clear d

% now we just need to set the fields for the output!
kinematicsAligned.time = kinTimeAligned;
kinematicsAligned.position = kinPosAligned;
kinematicsAligned.velocity = kinVelAligned;
kinematicsAligned.acceleration = kinAccAligned;

end

