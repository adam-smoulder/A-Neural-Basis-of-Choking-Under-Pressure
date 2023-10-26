function [dpcaInputMatrix] = formatDataForDPCA(data,labels,varargin)
% This function takes trial data (post stitching formatted) an dreorganizes
% it to be a good and tasty input for the dPCA toolbox.
%
% Inputs:
% - data: [ntrials x 1] cell vector with the [ntime x nfactors] factor
% scores in each cell.
% - labels: [ntrials x nlabels] matrix with the labels for this data. The
% number of labels (i.e., direction, reward) is the number of variables you
% want to demix (time is implicitly assumed to be used here). I'm not
% intending to try this for > 2 labels. 
% 
% Optional Inputs:
% - (1) downsampleFactor: scalar of factor by which to downsample time.
% Useful if you're using ms resolution data, which normally makes tensors
% that are too big for MATLAB. Default = 1
%
% Outputs:
% - dpcaInputMatrix: [nfactors x <nlabels> x ntime x maxNTrialsPerCond]
% tensor...what a mess.
%
% Adam Smoulder, 2/5/21

if isempty(varargin)
    downsampleFactor = 1;
else
    downsampleFactor = varargin{1};
end

% Get some variables
[ntrials,nlabels] = size(labels);
[ntime,nfactors] = size(data{1});

% Adjust ntime by the downsampleFactor
timeIndsToUse = 1:downsampleFactor:ntime;
ntimeDS = length(timeIndsToUse);

 % Convert labels to indices
 nconds_byLabel = nan(nlabels,1);
 for n = 1:nlabels
    curConds = unique(labels(:,n));
    [~,labels(:,n)] = ismember(labels(:,n),curConds);
    nconds_byLabel(n) = length(curConds);
 end; clear n

 % Count data by label condition
 [groupedData,inds_byCond] = groupDataByLabel(data,labels);
 ntrials_byCond = cellfun(@(x) length(x), inds_byCond);
 maxNTrialsPerCond = max(ntrials_byCond(:));
 
 % Loop over factors to get timepoints for each condition
 dpcaInputMatrix = nan([nfactors nconds_byLabel' ntimeDS maxNTrialsPerCond]);
 for f = 1:nfactors
     if nlabels == 1
         for n1 = 1:nconds_byLabel(1)
             curData = groupedData{n1};
              temp = cellfun(@(x) x(timeIndsToUse,f), curData, 'uniformoutput',false);
              dpcaInputMatrix(f,n1,:,1:ntrials_byCond(n1)) = [temp{:}];
         end; clear n1
     elseif nlabels == 2
         for n1 = 1:nconds_byLabel(1)
             for n2 = 1:nconds_byLabel(2)
                 curData = groupedData{n1,n2};
                 temp = cellfun(@(x) x(timeIndsToUse,f), curData, 'uniformoutput',false);
                 dpcaInputMatrix(f,n1,n2,:,1:ntrials_byCond(n1,n2)) = [temp{:}];
             end; clear n2
         end; clear n1
     elseif nlabels == 3
         for n1 = 1:nconds_byLabel(1)
             for n2 = 1:nconds_byLabel(2)
                 for n3 = 1:nconds_byLabel(3)
                    curData = groupedData{n1,n2,n3};
                    temp = cellfun(@(x) x(timeIndsToUse,f), curData, 'uniformoutput',false);
                    dpcaInputMatrix(f,n1,n2,n3,:,1:ntrials_byCond(n1,n2,n3)) = [temp{:}];
                 end; clear n3
             end; clear n2
         end; clear n1
     else
         error('This cannot handle > 3 labels yet...')
     end
 end; clear f
 
 
end

