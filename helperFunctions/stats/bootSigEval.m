function varargout = bootSigEval(data1,data2,alpha)
% bootSigEval assumes data1 and data2 are two bootstrapped distributions of
% equal length and determines if the difference between the data' means is
% significantly different from 0 at confidence (1-alpha). This is done
% using the "percentile bootstrap", meaning literally looing at if at least
% (1-alpha)*100% of the distribution of (data1-data2) lies on either side
% of 0.
%
% Note: 
% Because this is assumed to come from bootstrapped data, the input data
% is assumed to represent means, meaning that the CI is read directly
% from them (not normalized by the length of the data or anything)
%
% Inputs:
% data1 - input data, assumed to be means from a bootstrapped distribution
% data2 - same as data 1. The evaluation is then on data1-data2
% alpha - the significance level to test with the percentile bootstrap 
%
% Outputs:
% sig - 1 if significant, 0 if not
% p - the p value of the percentile bootstrap (so how much of the
% distribution of data1-data2 is past 0). This isn't really interpretable
% unless you have a huge amount of bootstraps.
%
% Adam Smoulder, 3/15/19

data = data1-data2;
fracAbove0 = sum(data > 0)/length(data);
fracBelow0 = sum(data <= 0)/length(data); % equality to 0 shouldn't matter

p = min([fracAbove0, fracBelow0]);
sig = p <= alpha;

switch nargout
    case 1
        varargout = sig;
    case 2
        varargout = [{sig},{p}];
    otherwise %...why
        throw(MException('MyComponent:inputError', 'Only 1 or 2 outputs allowed'));
end


end

