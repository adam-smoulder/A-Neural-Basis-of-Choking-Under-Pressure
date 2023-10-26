function [inputData, trialAssignments] = mapAssignments_Earl(inputData)
%%% This function takes an input dataset, inputData, which is assumed to be
%%% a structure of N trials in length with the internal structure
%%% 'Overview', along with a 'Parameters' structure and 'TrialData'. This
%%% is pretty much just meant for internal use on datasets from Aaron's
%%% lab.
%%% 
%%%
%%% inputs:
%%% - 'data':  see above for requirements
%%% - 'fname': name of data file used; should be fine so long as it
%%% includes the date of experiment
%%%
%%% outputs:
%%% - 'inputData':  same as input, though now includes 'code' field, which is a
%%% row vector of length R (however many identifying features are
%%% specified). 'code' is added at the top level. 'codeDefinition' is also
%%% added to trial 1 at this level (purely for reference) with the
%%% definition of what each index of 'code' reflects. Also, the size will
%%% now be R-a, where "a" is the number of trials in the data that don't 
%%% have a mapping in the root table (it's typically < 5)
%%% - 'assignments': N-a x R matrix of all codes, where N is the number of
%%% trials
%%% - 'lostTrials': indices of "a" mentioned above
%%% 
%%% Adam Smoulder, 9/26/18 (edit: 10/17/18)
%%% 

targetOffset = [90 -430];
allDirLabels = inputData(1).Parameters.endTargets(:,1:2).*[1 -1]-targetOffset;

% set stuff up
reachTarget_byTrial = nan(length(inputData),2); % 2nd dim is for X/Y coor
rewardSize_byTrial = nan(length(inputData),1);
dirLabel_byTrial = nan(length(inputData),1);
success_byTrial = nan(length(inputData),1);
trialLabels = nan(length(inputData),1);
stateTableLength = nan(length(inputData),1);

% for each trial, extract relevant info
for i = 1:size(rewardSize_byTrial,1) 
        % for later, get the stateTableLength (odd trials have 35 states)
        stateTableLength(i) = length(inputData(i).Parameters.StateTable);
        
        % get trial numbers
        trialLabels(i) = str2double(inputData(i).Overview.trialNumber(6:end));
        
        % Find the states used in the choice task.
        states = inputData(i).Parameters.stateNames;
        stateMatrix = nan(length(states),1);
        stateMatrix = num2cell(stateMatrix);
        for k = 1:size(stateMatrix,1)
            stateMatrix{k,1} = states{1,k};
        end
        
        % Find the reward size and target location associated with each target.
        targetHoldIndex = find(strcmp('Target Hold',stateMatrix) == 1);
        rewardPassState = inputData(i).Parameters.StateTable(targetHoldIndex).windowPassState;
        rewardSize_byTrial(i,1) = inputData(i).Parameters.StateTable(rewardPassState).Interval.length;           % Reward size of Reach

        % Get the target locations.
        reachTarget_byTrial(i,:) = inputData(i).Parameters.StateTable(targetHoldIndex).StateTargets.location(1,1:2).*[1 -1]-targetOffset;       % Target x and y of Reach
        
        % Give a direction label
        dirLabel_byTrial(i) = find(ismember(allDirLabels,reachTarget_byTrial(i,:),'rows'));
        
        % get success case
        success_byTrial(i) = inputData(i).Overview.trialStatus;
end

% convert reward sizes to reward labels
rewardSizes = unique(rewardSize_byTrial);
[~,rewardLabel_byTrial] = ismember(rewardSize_byTrial,rewardSizes);

% create trial assignments and update data and such
trialAssignments = [dirLabel_byTrial rewardLabel_byTrial success_byTrial trialLabels];
codeDefinition = [{'direction'},{'rewards'},{'success?'},{'trial number'}]; % if more parameters are added, edit this too!
for i = 1:length(inputData)
    inputData(i).codes = trialAssignments(i,:);
end
inputData(1).codeDefinition = codeDefinition;

end

