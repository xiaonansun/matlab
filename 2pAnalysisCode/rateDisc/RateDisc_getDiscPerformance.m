function [allPerf, h1] =  RateDisc_getDiscPerformance(bhv, ind, binSize, showPlot)

if ~exist('ind','var') || isempty(ind)
    ind = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %use all active trials
end

if ~exist('showPlot','var') || isempty(showPlot)
    showPlot = true;
end

if ~exist('binSize','var') || isempty(binSize)
    binSize = 10;
end

if ~isfield(bhv, 'flipIdx')
    bhv.flipIdx = false(1, length(ind));
end

%%
rInd = (bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded); %right-choice trials
rInd = rInd(ind);
rInd(bhv.flipIdx(ind)) = ~rInd(bhv.flipIdx(ind)); %flip left and rightward choices

tCnt = 0;
eventCnt = zeros(2,sum(ind));
for iTrials = find(ind)
    tCnt = tCnt +1;
    [left, right] = Behavior_getStimEvent(bhv.StimType(iTrials), bhv.stimEvents{iTrials});
    
    if ~bhv.flipIdx(iTrials)
        eventCnt(1,tCnt) = length(right);
    else
        eventCnt(1,tCnt) = length(left);
    end
    eventCnt(2,tCnt) = length(left) + length(right);
end

[nTrials, distRatio, trialInd]  = histcounts(eventCnt(1,:) ./ eventCnt(2,:), binSize); %get ratio between rightward and absolute number of events
distRatio = distRatio + diff(distRatio(1:2))/2; distRatio(end) = [];
for iBins = 1:length(nTrials)
    rightChoice(iBins) = sum(rInd(trialInd == iBins)); %get number of rightward trials for each difficulty
end
allPerf.distRatio = distRatio;

[params, h1, ~, cFit] = Behavior_fitPalamedes(distRatio, rightChoice, nTrials, showPlot, false); %complete discrimination parameters
allPerf.cFit = cFit;

fNames = fieldnames(params);
for iFields = 1 : length(fNames)
    allPerf.(fNames{iFields}) = params.(fNames{iFields});
end

% discrimination only
if sum(rightChoice(2:end-1)) > 0
    [params, ~, ~, cFit] = Behavior_fitPalamedes(distRatio(2:end-1), rightChoice(2:end-1), nTrials(2:end-1)); %complete discrimination parameters
    allPerf.disc_cFit = cFit;
    if showPlot
        hold(h1.fit.Parent, 'on');
        h2 = plot(h1.fit.Parent,cFit(1,:),cFit(2,:),'color','r', 'linewidth', 4);
        uistack(h2,'bottom');
        hold(h1.fit.Parent, 'off');
    end
    
    fNames = fieldnames(params);
    for iFields = 1 : length(fNames)
        allPerf.(['disc_' fNames{iFields}]) = params.(fNames{iFields});
    end
end

allPerf.disc = rightChoice ./ nTrials;
[upper, lower] = Behavior_wilsonError(rightChoice, nTrials);
allPerf.discUpper = upper;
allPerf.discLower = lower;