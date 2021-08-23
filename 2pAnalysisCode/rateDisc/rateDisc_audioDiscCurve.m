function [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh, stats] = rateDisc_audioDiscCurve(bhv, cInd, distBins, discOnly, fixBias, returnCIs)
% short code to compute auditory discrimination performance for a subset of trials in the bhv structure.


if ~exist('distBins','var') || isempty(distBins)
    distBins = 10;
end
if ~exist('discOnly','var') || isempty(discOnly)
    discOnly = false;
end
% optional input to fix the bias term
if ~exist('fixBias','var') || isempty(fixBias)
    fixBias = false; %this fixes the bias (x-axis shift) centered between left and right
end
% optional input to return boot-strapped statistic for model estimates
if ~exist('returnCIs','var') || isempty(returnCIs)
    returnCIs = false;
end

ind = cInd & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.StimType == 2; %only use active audio trials
rInd = (bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded); %right-choice trials
rInd = rInd(ind);

tCnt = 0;
eventCnt = zeros(2,sum(ind));
for iTrials = find(ind)
    tCnt = tCnt +1;
    [left, right] = Behavior_getStimEvent(bhv.StimType(iTrials), bhv.stimEvents{iTrials});
    eventCnt(1,tCnt) = length(right);
    eventCnt(2,tCnt) = length(left) + length(right);
end

[nTrials, distRatio, trialInd]  = histcounts(eventCnt(1,:) ./ eventCnt(2,:), distBins); %get ratio between rightward and absolute number of events
distRatio = distRatio + diff(distRatio(1:2))/2; distRatio(end) = [];
for iBins = 1:length(nTrials)
    rightChoice(iBins) = sum(rInd(trialInd == iBins)); %get number of rightward trials for each difficulty
end

pChoseHigh = rightChoice./nTrials;
if ~any(rightChoice == 0)
    if discOnly
        [params, ~, stats, cFit] = Behavior_fitPalamedes(distRatio(2:end-1), rightChoice(2:end-1), nTrials(2:end-1), false, true, fixBias, returnCIs); %complete discrimination parameters
    else
        [params, ~, stats, cFit] = Behavior_fitPalamedes(distRatio, rightChoice, nTrials, false, fixBias, returnCIs); %complete discrimination parameters
    end
    
    %Get 95% Wilson binomial CIs for data
    z = 1.96;
    dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
    dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
else
    params = [];
    cFit = [];
    dataUpper = []; dataLower = [];
end