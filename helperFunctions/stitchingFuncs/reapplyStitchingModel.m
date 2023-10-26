function [model,lEstO,missingStuffEsts] = reapplyStitchingModel(data,model)
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
% Adam Smoulder, 3/4/19

disp('Beginning reapplying stitching model!')

% Infer latents
lEst = inferFALatents(model, data);

% Put latents into an orthonormalized coordinate system
[model, lEstO] = orthonormalizeFAMdl(model, lEst);

% Predict missing data
missingStuffEsts = data;
nDataPts = size(data,1);

for i = 1:nDataPts
    curData = data(i,:);
    missingInds = isnan(curData);
    missingPred = model.d(missingInds) + model.C(missingInds,:)*lEst(i,:)';
    missingStuffEsts(i,missingInds) = missingPred;
end

disp('Completed reapplying stitching model!')

end

