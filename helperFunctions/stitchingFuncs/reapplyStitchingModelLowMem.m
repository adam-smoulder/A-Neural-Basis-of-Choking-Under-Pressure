function [model,lEstO] = reapplyStitchingModelLowMem(data,model,nsplit)
% Same as reapplyStitchingModel.m but doesn't output the estimates of
% missing values and runs nsplit number of segments of time bins
% individually (e.g. if you had 10^7 time bins and used nsplits = 100, it
% would divide the data into sets of 10^5 time bins to reapply the model
% to). For even more memory saving, change the input data to a single
% instead of a double
%
%
% This function is simply a c&p of the end of faCompleteDataMatrix used to
% apply a stitched FA model to further data. This is useful for instances
% such finding the latents for data that has been realigned or projecting
% new data into the same model.
%
% Inputs: 
% data - the data to input - in the same format as if it were going into
% stitching. time x neurons
% model - the existing stitching model to reapply
%
% Outputs:
% model - same as the input
% lEstO - orthonormalized latents from the data (time x latents)
% missingStuffEsts - the original data matrix with the NaNs filled in
%
%
% Adam Smoulder, 6/21/19

disp('Beginning reapplying stitching model (with low mem)!')

% First, we split the data into nsplit separate timebin sets to run individually
timeSplits = ceil(linspace(0,size(data,1),nsplit+1)); % start at 0 for indexing more easily in the loop
lEstO = nan(size(data,1),size(model.C,2));

% Then we sequentially reapply the model to each!
for s = 1:nsplit   
    curInds = (timeSplits(s)+1):timeSplits(s+1);
    curData = data(curInds,:);
    
    % Infer latents
    lEst = inferFALatents(model, curData);
    
    % Put latents into an orthonormalized coordinate system
    [~, lEstOTemp] = orthonormalizeFAMdl(model, lEst);
    lEstO(curInds,:) = lEstOTemp;
    
    disp(['Reapplied model for segment ' num2str(s) ' of ' num2str(nsplit)])
end; clear s

disp('Completed reapplying stitching model!')

end

