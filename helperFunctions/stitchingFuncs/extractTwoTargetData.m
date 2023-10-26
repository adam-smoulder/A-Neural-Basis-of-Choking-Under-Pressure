function [data,selectedTrialAssignments,allTrialAssignmentOptions] = extractTwoTargetData(data)
% This function takes in data that are assumed to be two-target trials (all
% others are removed), extracts some of their relevant info for both 
% options (i.e. tubeChoicesLabels, rewardChoicesLabels) and decided
% quantities (i.e. tubeLabels, rewardLabels). Punishments are not processed
% here but are used in the case of false starts to assume intended tube
% selection.
%
% Input/Output: 
% - data: [N x 1] structure, where N is the number of two-target trials.
% This function adds the field "TrialParams" onto data, which contains the
% reward/tube/punish choices, the selected reward/tube/punish, and what
% option was selected
%
% Outputs:
% - selectedTrialAssignments: [N x 5] matrix containing [tube reward punish success trial]
% labels for each trial. This is the same format as the mapAssignments
% method's output!
% - allTrialAssignmentOptions: Structure with the choices for each trial. I
% made a structure out of laziness, which I may regret later. Includes
% fields tubeChoices, rewardChoices, punishChoices, directionChoices, and
% selections
%
% Adam Smoulder, 4/18/19 (edit 4/22/19)


% find first multitarget trial and only keep that data
firstMultitarg = 0;
counter = 1;
while firstMultitarg == 0
    if(length(data(counter).Parameters.StateTable)>= 15)
        firstMultitarg = counter;
    end
    counter = counter+1;
end

data = data(firstMultitarg:end);


% set stuff up
reachTarget_choice1 = nan(length(data),2); % 2nd dim is for X/Y coor
reachTarget_choice2 = nan(length(data),2);
rewardSize_byTrial = nan(length(data),2); % 2nd dim is for choice 1 vs 2
punishState_byTrial = nan(length(data),2);
tubeChoicesLabels = nan(length(data),2);
success = nan(length(data),1);
choiceSelections = nan(length(data),1);
trialLabels = nan(length(data),1);
stateTableLength = nan(length(data),1);

% for each trial, extract relevant info
for i = 1:size(rewardSize_byTrial,1) 
        % for later, get the stateTableLength (odd trials have 35 states)
        stateTableLength(i) = length(data(i).Parameters.StateTable);
        
        % get trial numbers
        trialLabels(i) = str2double(data(i).Overview.trialNumber(6:end));
        
        % Find the states used in the choice task.
        states = data(i).Parameters.stateNames;
        stateMatrix = nan(length(states),1);
        stateMatrix = num2cell(stateMatrix);
        for k = 1:size(stateMatrix,1)
            stateMatrix{k,1} = states{1,k};
        end
        
        % Find the reward size and target location associated with each target.
        targetHoldIndex = find(strcmp('Target Hold',stateMatrix) == 1);
        targetHoldIndex1 = find(strcmp('Target Hold1',stateMatrix) == 1);
        
        rewardPassState = data(i).Parameters.StateTable(targetHoldIndex).windowPassState;
        rewardPassState1 = data(i).Parameters.StateTable(targetHoldIndex1).windowPassState;
        
        rewardSize_byTrial(i,1) = data(i).Parameters.StateTable(rewardPassState).Interval.length;           % Reward size of Reach
        rewardSize_byTrial(i,2) = data(i).Parameters.StateTable(rewardPassState1).Interval.length;          % Reward size of Reach 1

        punishState_byTrial(i,1) = data(i).Parameters.StateTable(targetHoldIndex).windowFailState;        % punish state of Reach
        punishState_byTrial(i,2) = data(i).Parameters.StateTable(targetHoldIndex1).windowFailState;       % punish state of Reach 1

        % Get the target locations.
        reachTarget_choice1(i,:) = data(i).Parameters.StateTable(targetHoldIndex).StateTargets.location(1,1:2).*[1 -1]-[100 -510];       % Target x and y of Reach
        reachTarget_choice2(i,:) = data(i).Parameters.StateTable(targetHoldIndex1).StateTargets.location(1,1:2).*[1 -1]-[100 -510];
        
        % Get the tubes presented (left to right)
        tubeTraj = data(i).Parameters.TrialTubeParameters.trajectory(:,1:2).*[1 -1]-[100 -510];
        tubeSigns = sign(tubeTraj([25 75],:)); % [x y ; x y] for tube midpoints
        tubeYSigns(1) = tubeSigns((tubeSigns(:,1)==-1),2); % left tube
        tubeYSigns(2) = tubeSigns((tubeSigns(:,1)==1),2);  % right tube
        switch tubeYSigns(1)*10+tubeYSigns(2)
            case -11 % w shaped (down, down)
                tubeChoicesLabels(i,:) = [6 1];
            case -9 % negative sine curve (down, up)
                tubeChoicesLabels(i,:) = [6 2];
            case 9  % positive sine curve (up, down)
                tubeChoicesLabels(i,:) = [5 1];
            case 11 % m shaped (up, up)
                tubeChoicesLabels(i,:) = [5 2];
            otherwise
                disp('FUCK we shouldnt be seeing this')
                tubeYSigns
        end
        
        % Currently, the tubes are ordered such that 5/6 is "choice 1" and
        % 1/2 is "choice 2". To correct this, we can identify which reach
        % target is choice 1 vs. 2 from the reachTarget variables. If the
        % choice 1 reachTarget is leftward (-), then tube 5/6 should be
        % associated with choice 1, and all is good. If the choice 1
        % reachTarget is rightward (+), we need to flip the tube order.
        if reachTarget_choice1(i,1)==100
            tubeChoicesLabels(i,:) = fliplr(tubeChoicesLabels(i,:));
        end
        
        
        
        % identify which option was selected (0 if it doesn't make it to
        % selection)
        sts = data(i).TrialData.stateTransitions;
        if sum(sts(1,:)==8)==1 % first choice
            choiceSelections(i) = 1;
        elseif sum(sts(1,:)==9)==1 % 2nd choice
            choiceSelections(i) = 2;
        else % either we didn't make it to the reach or it just started. Let's look at the velocity to check this
            if sum(sts(1,:)==7) ~= 0 % if we can get to the reach we can see which way he moves
                reachTime = double(sts(2,find(sts(1,:)==7,1,'first')));
                kinStruct = data(i).TrialData.HandKinematics.markerKinematics;
                position = kinStruct.position(:,1:2).*[1 -1]-[100 -510];
                reachInd = find((kinStruct.time > reachTime),1,'first');
                startPos = position(reachInd-40,:);
                pos150ms = position(reachInd+20,:);
                x_displacement = pos150ms(1)-startPos(1);
                if x_displacement > 2 % he chose right
                    if reachTarget_choice1(i,1) == 100 % choice 1 is right, choice 2 is left
                        choiceSelections(i) = 1;
                    else % choice 1 is left, choice 2 is right
                        choiceSelections(i) = 2;
                    end
                elseif x_displacement < -2 % he chose left
                    if reachTarget_choice1(i,1) == 100 % choice 1 is right, choice 2 is left
                        choiceSelections(i) = 2;
                    else % choice 1 is left, choice 2 is right
                        choiceSelections(i) = 1;
                    end
                else
                    choiceSelections(i) = 0; % not enough displacement, so we don't know :(
                end
            else % if no reach, we can't say which tube
                choiceSelections(i) = 0; 
            end
        end
        
        % get success case
        success(i) = data(i).Overview.trialStatus;
end
directionChoicesLabels = 2-(round([reachTarget_choice1(:,1) reachTarget_choice2(:,1)])/100);


% process extracted data:
% % correct tube directions based on the target labels THIS IS DONE EARLIER NOW!!!
% for i = 1:length(data)
%     if directionChoicesLabels(i,1)==1 % the right tube is on the first slot...
%         tubeChoicesLabels(i,:) = fliplr(tubeChoicesLabels(i,:));
%     end
% end; clear i

% map rewards to easy values (1-4)
rewardMagnitudes = unique(rewardSize_byTrial);
assert(length(rewardMagnitudes)==4)
rewardConverter(rewardMagnitudes)=[1 2 3 4];
rewardChoicesLabels = rewardConverter(rewardSize_byTrial);

% map punishments to easy values
oddTrials = find(stateTableLength~=25); % for some odd reason, some trials have more states...?
tempPunishChoicesLabels = punishState_byTrial-15;
punishChoicesLabels = tempPunishChoicesLabels;
punishChoicesLabels(oddTrials',:) = round((tempPunishChoicesLabels(oddTrials',:)+0.1)/2); % these need updated because there's extra states after 15


% Identify labels based on the selection made
% first, guess which decision would've been made based on optimal r/p if
% none was made
rewardLabels = nan(size(rewardChoicesLabels(:,1)));
tubeLabels = nan(size(rewardChoicesLabels(:,1)));
directionLabels = nan(size(rewardChoicesLabels(:,1)));
punishLabels = nan(size(rewardChoicesLabels(:,1)));

noSelectionTrials = (choiceSelections==0);
selectionTrialInds = find(~noSelectionTrials);
noSelectionTrialInds = find(noSelectionTrials);
for i = 1:length(selectionTrialInds)
    rewardLabels(selectionTrialInds(i)) = rewardChoicesLabels(selectionTrialInds(i),choiceSelections(selectionTrialInds(i)));
    tubeLabels(selectionTrialInds(i)) = tubeChoicesLabels(selectionTrialInds(i),choiceSelections(selectionTrialInds(i)));
    directionLabels(selectionTrialInds(i)) = directionChoicesLabels(selectionTrialInds(i),choiceSelections(selectionTrialInds(i)));
    punishLabels(selectionTrialInds(i)) = punishChoicesLabels(selectionTrialInds(i),choiceSelections(selectionTrialInds(i)));
end; clear i
for i = 1:length(noSelectionTrialInds) % if false start (no selection made)
    if diff(rewardChoicesLabels(noSelectionTrialInds(i),:),[],2) == 0 % if same reward size
        if diff(punishState_byTrial(noSelectionTrialInds(i),:),[],2) == 0 % if same punish size too, it doesn't matter
            rewardLabels(noSelectionTrialInds(i)) = rewardChoicesLabels(noSelectionTrialInds(i),1);
            tubeLabels(noSelectionTrialInds(i)) = 0; % no tube decided;
            directionLabels(noSelectionTrialInds(i)) = 0;
            punishLabels(noSelectionTrialInds(i)) = punishChoicesLabels(noSelectionTrialInds(i),1);
        else % pick lower punishment
            [~,sel] = min(punishState_byTrial(noSelectionTrialInds(i),:));
            rewardLabels(noSelectionTrialInds(i)) = rewardChoicesLabels(noSelectionTrialInds(i),sel);
            tubeLabels(noSelectionTrialInds(i)) = tubeChoicesLabels(noSelectionTrialInds(i),sel);
            directionLabels(noSelectionTrialInds(i)) = directionChoicesLabels(noSelectionTrialInds(i),sel);
            punishLabels(noSelectionTrialInds(i)) = punishChoicesLabels(noSelectionTrialInds(i),sel);
        end
    else % pick greater reward
        [~,sel] = max(rewardChoicesLabels(noSelectionTrialInds(i),:));
        rewardLabels(noSelectionTrialInds(i)) = rewardChoicesLabels(noSelectionTrialInds(i),sel);
        tubeLabels(noSelectionTrialInds(i)) = tubeChoicesLabels(noSelectionTrialInds(i),sel);
        directionLabels(noSelectionTrialInds(i)) = directionChoicesLabels(noSelectionTrialInds(i),sel);
        punishLabels(noSelectionTrialInds(i)) = punishChoicesLabels(noSelectionTrialInds(i),sel);
    end
end; clear i


% reassign relevant info to data and outputs
selectedTrialAssignments = [tubeLabels rewardLabels punishLabels success trialLabels];
allTrialAssignmentOptions.tubeChoices = tubeChoicesLabels;
allTrialAssignmentOptions.rewardChoices = rewardChoicesLabels;
allTrialAssignmentOptions.punishChoices = punishChoicesLabels;
allTrialAssignmentOptions.directionChoices = directionChoicesLabels;
allTrialAssignmentOptions.selections = choiceSelections;
for i = 1:length(data)
    data(i).TrialParams.selections = choiceSelections(i);
    data(i).TrialParams.rewardChoices = rewardChoicesLabels(i,:);
    data(i).TrialParams.reward = rewardLabels(i);
    data(i).TrialParams.punishChoices = punishChoicesLabels(i,:);
    data(i).TrialParams.punish = punishLabels(i);
    data(i).TrialParams.tubeChoices = tubeChoicesLabels(i,:);
    data(i).TrialParams.tube = tubeLabels(i);
    data(i).TrialParams.directionChoices = directionChoicesLabels(i,:);
    data(i).TrialParams.direction = directionLabels(i);
    data(i).TrialParams.trialNumber = trialLabels(i);
    data(i).TrialParams.success = success(i);
    data(i).codes = selectedTrialAssignments(i,:);
end; clear i
data(1).codeDefinition = [{'tubes'},{'rewards'},{'punishments'},{'success?'},{'trial number'}]; % if more parameters are added, edit this too!

end

