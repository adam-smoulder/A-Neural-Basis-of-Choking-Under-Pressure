function [stitchingInput] = stitchingPrep_appendLag(stitchingInput,nlags)
%UNTITLED31 Summary of this function goes here
%   Detailed explanation goes here
X = stitchingInput.obs;
lastTimeBins = stitchingInput.lastTimeBins;
ndays = length(lastTimeBins);

megaX = [];

for j = 1:nlags
    XplanLag = nan(size(X));
    for i = 1:ndays % we have to do it by day or else the nans will get us
        if i == 1
            XplanLag(1:lastTimeBins(1),:) = [zeros([j,size(X,2)]) ; X(1:lastTimeBins(1)-j,:)];
        else
            XplanLag((lastTimeBins(i-1)+1):lastTimeBins(i),:) = [zeros([j,size(X,2)]) ; X((lastTimeBins(i-1)+1):(lastTimeBins(i)-j),:)];
        end;
    end
    megaX = [XplanLag megaX];
end;


stitchingInput.obs = [megaX X];
disp(["Appended " num2str(nlags) " lag(s) to data"])
