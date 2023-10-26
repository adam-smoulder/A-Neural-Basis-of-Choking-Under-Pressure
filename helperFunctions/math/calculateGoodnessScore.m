function [goodnessScore,distFromBoundary,niterToComplete] = calculateGoodnessScore(input,labels)
% See calculateGoodnessScore_trainTest; this is just a breakout of it where
% you train and test with the same data

[goodnessScore,distFromBoundary,niterToComplete] = calculateGoodnessScore_trainTest(input,labels,input,labels);

end

