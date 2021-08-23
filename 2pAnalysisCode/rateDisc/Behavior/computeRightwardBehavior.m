function [perf,error,tCount] = computeRightwardBehavior(bhv,maxCount,varargin)
% Function to compute behavioral performance, given a set of constraints to
% select a subset of trials from a larger bhv structure. Returns
% rightward choices as a fraction between 0 and 1 and an error term which 
% is the 95% confidence interval assuming a binomial distribution. 
% Also returns the absolute amount of included trials in each condition.
%
% Usage: [perf,error,tCount] = selectBehavior(bhv,varargin)
% 
% 'bhv' is a behavioral structure that is created from Bpod and containts 
% two fields 'Rewarded' and 'Punished' that are used to compute performance.
% maxCount is the maximal amount of trials that is used to compute 
% performance. If more trials are available, only the first 'maxCount'
% trials are being used.
%
% Every first field in 'varargin' is string that points to a field in the bhv
% structure and has to be followed by either a number or string that can be
% used to select trials from that field. 
% For example, computeBehavior(bhv,'TargStim',15) will select all trials that
% contain the value 15 for the field TargStim.
% Alternatively, varargin can be a single logical index with the same
% length as trials in the bhv structure. In this case the index defines the
% used trials directly.

%% set up constraints and find subset of trials
if length(varargin) > 1
    tCheck = [];
    for x = 1:2:length(varargin);
        tCheck = [tCheck 'bhv.' varargin{x} ' == ' num2str(varargin{x+1})];
        if x ~= length(varargin)-1
            tCheck = [tCheck ' & '];
        end
    end
    eval(['ind = ' tCheck ';']);
else
    ind = varargin{1};
end

ind = find(ind);

if isempty(maxCount)
    maxCount = 500;
end
if length(ind) > maxCount
   ind = ind(1:maxCount);
end

tCount = sum(bhv.Rewarded(ind)+bhv.Punished(ind)); %trialcount
perf = sum(bhv.Rewarded(ind))/tCount; %peformance for given distractor
error = binoinv([0.025 0.975], tCount, perf) ./ tCount; %convidence intervals for given distractor
       
