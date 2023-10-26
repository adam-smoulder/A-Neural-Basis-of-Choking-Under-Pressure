function [outlierLabels] = identifyOutliersInNeuralData(data,dayLabels,displayIterations)
% This function identifies outliers in data combined across days in the
% following manner:
% - Receive the data across days as input (data, ntrials x nneurons)
% - For each day (gleaned from dayLabels input):
%   1. Extract the data for the current day
%   2. Do PCA and get out PC1
%   3. Remove outliers with high criterion (7.5*MAD, moving 100 trial
%   window)
%   4. Repeat steps 2-3 until no outliers remain
% Any trial that exceeds the 7.5*MAD criterion are labeled as bad.
%
% Inputs:
% - data: [ntrials x nneurons] firing rate data, presumably the firing rate
% for some bin taken from each trial. This is combined across days, with
% nans being used for indices when a neuron is not present for a given
% trial.
% - dayLabels: [ntrials x 1] label indicating which day the trial is from.
% - displayIterations: boolean, if true shows figure for each iteration
%
% Outputs:
% - outlierLabels: [ntrials x 1] boolean of if the trial is an outlier or
% not.
%
% Adam Smoulder, 2/16/23

ntrials = size(data,1);
outlierLabels = false(ntrials,1);
days = unique(dayLabels); ndays = length(days);
for a = 1:ndays
    % Get the current day's data
    curInds = find(dayLabels==days(a));
    curData = data(curInds,:);
    curData(:,isnan(sum(curData))) = [];
    curntrials = size(curData,1);
    
    % For this day, check for outliers with a 100 trial moving median and
    % 6*median absolute deviation criterion against the norm of the top 3
    % PCs
    curValidInds = (1:curntrials)';
    noutliers_rep = inf; % just to start the loop
    iter = 1;
    while noutliers_rep ~= 0
        [~,zOutTest] = pca(curData(curValidInds,:),'numcomponents',3);
        zOutTest = sqrt(sum(zOutTest.^2,2));
        outTestMovMed = movmedian(zOutTest,101,'endpoints','discard');
        outTestMovMAD = movmad(zOutTest,101,'endpoints','discard');
        outTestMovThresh = outTestMovMed+6*outTestMovMAD; % median equiv of mean + 4 std
        outTestMovThresh = [outTestMovThresh(1)*ones(50,1) ; outTestMovThresh ; outTestMovThresh(end)*ones(50,1)];
        repOutliers = find(zOutTest > outTestMovThresh);
        curValidInds(repOutliers) = [];
        noutliers_rep = length(repOutliers);
        
        % Show figures
        if displayIterations
            figure; plot(zOutTest,'o'); hold on ; plot(repOutliers,zOutTest(repOutliers),'*','linewidth',2); plot(outTestMovThresh,'k-')
            title(['Day ' num2str(a) ' iter ' num2str(iter)]); ylabel('L2 norm of PC 1-3'); xlabel('Trial #')
            pause(0.25)
        end
        
        iter = iter+1;
    end
    curOutlierLabels = ~ismember((1:curntrials)',curValidInds);
    outlierLabels(curInds(curOutlierLabels)) = true;
end; clear a

end

