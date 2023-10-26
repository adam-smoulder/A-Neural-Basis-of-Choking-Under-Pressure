function [outputData, assignments, lostTrialIndices] = mapAssignments(data, fname)
%%% createAssignments takes an input dataset, data, which is assumed to be
%%% a structure of N trials in length with the internal structure
%%% 'Overview'. Also uses the trialDescriptionData.txt table 
%%% 
%%% 'Overview' should contain a field, 'trialNumber', which contains a
%%% string "Trial##". The ## will be used for comparison with the table of
%%% info from trialDescriptionData.txt.
%%%
%%% inputs:
%%% - 'data':  see above for requirements
%%% - 'fname': name of data file used; should be fine so long as it
%%% includes the date of experiment
%%%
%%% outputs:
%%% - 'data':  same as input, though now includes 'code' field, which is a
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


% define parameters we care about
% tubes = 1:8;
% rewards = 1:4;
% punishments = 1:4;
% success = [0 1];
codeDefinition = [{'tubes'},{'rewards'},{'punishments'},{'success?'},{'trial number'}]; % if more parameters are added, edit this too!

% get some stuff we'll need...
trialInfo = readtable('trialDescriptionData.txt');
dateIndexInFname = strfind(fname,'20');                                        % should work till we collect data from 2100+ or change naming conventions 
dateOfExperiment = fname(dateIndexInFname(1):(dateIndexInFname(1)+7));         % we choose the 1st element incase 20 comes up again in fname
dateOfExperimentIndices = find(trialInfo{:,1}==str2double(dateOfExperiment));
trialInfoForData = trialInfo(dateOfExperimentIndices,:);                       % isolate trials for the given date

% get trial numbers from data and from the table
trialNumsFromData = []; 
for i = 1:length(data)
    trialNumsFromData = ...
        [trialNumsFromData ; str2double(data(i).Overview.trialNumber(6:end))]; 
end; clear i;
trialNumsFromTable = trialInfo{dateOfExperimentIndices,3};

% if it's in the data but not in the table, remove it
[~, lostTrialIndices] = setdiff(trialNumsFromData, trialNumsFromTable);
data(lostTrialIndices) = [];
trialNumsFromData(lostTrialIndices) = [];
ntrials = length(data);

outputData = [];
assignments = [];
for i = 1:length(trialNumsFromTable) % get info for each table row
    % find table trial num that aligns with data trial num
    dataIndex = find(trialNumsFromData == trialNumsFromTable(i));
    if isempty(dataIndex)
        continue % if it's in the table but not the data, skip it
    elseif length(dataIndex) > 1
        dataIndex = dataIndex(1); % pick the first... this shouldn't happen
        disp('Two trial indices for one data index...check that dates were separated properly')
    end
    
    % update the output (new) data and assignments
    trialTube = tubeRemapping(trialInfoForData{i,7}); % this part is bulky but done for clarity
    trialReward = trialInfoForData{i,8};
    trialPunishment = trialInfoForData{i,9};
    trialSuccess = trialInfoForData{i,5};
    trialNumber = trialInfoForData{i,3}; 
    assert(trialNumber==trialNumsFromTable(i)); % just for confirmation...
    
    trialAssignments = [trialTube trialReward trialPunishment trialSuccess trialNumber];
    data(dataIndex).codes = trialAssignments;
    outputData = [outputData ; data(dataIndex)];     % append trial to the output
    assignments = [assignments ; trialAssignments];  % append code for assignments
end; clear i

outputData(1).codeDefinition = codeDefinition;

end

