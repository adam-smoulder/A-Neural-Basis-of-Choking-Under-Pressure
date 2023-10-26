function [w,projData,eigVals] = orthogonalLDA(data,labels)
% The goal of this function is to give a projection matrix that both
% maximizes linear separation of data for the provided class labels and
% also enforces the projection dimensions to be orthonormal. Standard LDA
% does not enforce orthogonality; we do so here by simply finding
% discriminant dimensions one-by-one and subtracting the up-projected
% result of each LDA from our data before the next LDA is performed.
%
% ASSUMES NUMBER OF CLASSES (K) IS LESS THAN NUMBER OF OBSERVATIONS (N)
%
% Inputs:
% - data: N x D matrix of data with N observations and D dimensions
% - labels: N x 1 vector of labels for the data with K unique classes
%
% Outputs:
% - w:  D x K-1 orthonormal projection matrix, ordered by discriminability
% - projData:  the data projected onto the orthogonal LDA plane, (data*w)
% - eigVals:  Eigenvalues from the generalized eigenproblem solved for LDA;
% gives the rough "importance" of each of the output dimensions in
% discrimination of the given classes
%
% It should be noted that, unlike normal LDA, this method can be
% technically continued into perpetuity beyond K-1 times, as the solution
% is iterative (not closed form). No clue what you get if you keep
% going...I doubt anything of substance. Feel free to try.
%
% Adam Smoulder, 9/24/19

data = data-mean(data); % mean center data
classes = unique(labels);
K = length(classes);
D = size(data,2);

w = nan(D,K-1);
eigVals = nan(K-1,1);
curData = data; % for us to edit as we go
for i = 1:(K-1)
    discMdl = fitcdiscr(curData,labels);
    sB = discMdl.BetweenSigma;
    sW = discMdl.Sigma;
    [curEigVecs, curEigVls] = eig(sB, sW); % LDA dims and their relative importance
    
    % We only want the best dim for discrimination in each step
    [maxEigVal, bestEigVecInd] = max(diag(curEigVls));
    eigVals(i) = maxEigVal;
    curW = curEigVecs(:,bestEigVecInd);
    
    % This vector can have arbitrary values in the directions that have
    % already been removed; we subtract the projections on these directions
    % and normalize to ensure the output w is orthonormal
    if i > 1
        for j = 1:(i-1)
            oldW = w(:,j);
            curW = curW-(oldW'*curW)*oldW; % subtract out the projection in the irrelevant (collapsed) dimension
        end; clear j
    end
    curW = curW/norm(curW); % it's a vector, so norm works fine
    
    % Append this to our output matrix and subtract the up projection from
    % the data
    curData = curData - curData*(curW*curW');
    w(:,i) = curW;
end; clear i

% Finally, project our data onto this space
projData = data*w;
end

