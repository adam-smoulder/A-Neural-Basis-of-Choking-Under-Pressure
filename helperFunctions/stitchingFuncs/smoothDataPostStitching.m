function [lEstOSmooth,masterInfo] = smoothDataPostStitching(lEstO,masterInfo,kernel,doOverlap,varargin)
% This function is made to smooth (or rebin) data after stitching is done.
% This typically assumes the data uses small time bins (e.g. 5ms), though
% it can technically be done with larger ones. 
%
% Inputs:
% - lEstO: [ntimebins x nfactors] Matrix of orthonormalized factor scores.
% - masterInfo: [ntimebins x nlabels] Matrix with many different labels. If
% no overlap, certain values will be removed. We use the first two columns
% to apply the kernel on each day-trial. 4th column (time) used as well if
% no overlap is selected.
% - kernel: [nsmooth x 1] Kernel to be used to smooth the data. Will be
% normalized to sum to 1.
% - doOverlap: boolean, do we leave overlapping bins (true) or do we remove
% them (false)?
%
% Optional Inputs:
% - doFast: boolean, should we do this the fast and dirty way (not trial by
% trial, but just on the full matrix at once (true)? Default is false. This 
% assumes that doOverlap is false. The upside is that this is way faster
% than convolving trial-by-trial carefully. The downside is that you then
% will have edge effects. If you're not near the edges or have a small
% kernel, this should be safe and much less painful than the normal method.
%
% Outputs:
% - lEstO: same as input but smoothed
% - masterInfo: same as input but with bins removed if no overlap
%
% Some quick tips:
% - If no overlap, kernel length will determine the new bin size.
% - To do "rebinning", just make the kernel all ones for the length desired
% to rebin and set doOverlap to false.
% 
% Adam Smoulder, 4/15/20

doFast = false;
if length(varargin)~=0
    doFast = varargin{1};
end

% First, assert kernel is a column vector and normalize it
assert(sum(size(kernel))-1 == length(kernel)) % assert 1D
if length(kernel)==size(kernel,2) % if it's a row
    kernel = kernel';
end
kernel = kernel/sum(kernel);

% Apply kernel to all data
lEstOSmooth = nan(size(lEstO));
dayTrialLabels = masterInfo(:,1)*10000+masterInfo(:,2);
dayTrials = unique(dayTrialLabels);
indsToKeep = zeros(size(lEstO,1),1);

if ~doFast % do it trial-by-trial
    for i = 1:length(dayTrials) % for each trial
        curInds = dayTrialLabels==dayTrials(i);
        
        % Actually do the smooth
        for j = 1:size(lEstO,2) % for each variable
            lEstOSmooth(curInds,j) = conv(lEstO(curInds,j),kernel,'same'); % zero-pad the edges.
        end; clear j
%         if i == 1, figure; hold on; plot(lEstO(curInds,1)); plot(lEstOSmooth(curInds,1)); end
        
        % If nonoverlapping bins should be removed, we do it.
        curKeepInds = zeros(sum(curInds),1);
        if ~doOverlap
            % Find the bin closest to 0; we keep this one first
            curTime = masterInfo(curInds,4);
            [~,closestZeroInd] = min(abs(curTime));
            
            % For inds to keep, we take steps based on the size of the kernel
            % from the closestZeroInd
            posIndsToKeep = closestZeroInd:length(kernel):length(curInds);
            negIndsToKeep = fliplr(closestZeroInd:length(kernel):1);
            curIndsToKeep = [negIndsToKeep posIndsToKeep(2:end)]; % make sure not to double count the closestZeroInd
            
            % Mark which indices to keep
            curKeepInds(curIndsToKeep) = 1;
        else % keep all of em
            curKeepInds = 1+curKeepInds;
        end
        indsToKeep(curInds,:) = curKeepInds;
        
        % display progress sometimes
        if mod(i,50)==0
            disp(['Smoothed trial ' num2str(i) ' of ' num2str(length(dayTrials))])
        end
    end; clear i
    
else % if doFast, we go allll at once
    disp('Applying kernel to all data at once (will have trial edge overlaps)')
    indsToKeep = ones(size(lEstO,1),1);
    for j = 1:size(lEstO,2) % for each variable
        lEstOSmooth(:,j) = conv(lEstO(:,j),kernel,'same');
    end; clear j
end

% Keep only the indices marked to be kept
lEstOSmooth = lEstOSmooth(indsToKeep==1,:);
masterInfo = masterInfo(indsToKeep==1,:); 

disp('Done smoothing stitched output')

figure; hold on; plot(lEstO(1:500,1)); plot(lEstOSmooth(1:500,1));
end

