function allData = rateDisc_taskEpisodesHistoryOpto(bhv,cInd,dateRange)
%short function to compute detection performance with optogenetic
%inactivations during different task episodes and cortical locations.
%this takes into account trial history by analyzing trials separately,
%depending on whether the previous trial was rewarded or not.

typeOrder = [5 1 4 2 3]; %order different task episodes chronologically. This assumes Handle(5), Stimulus(1), LateStimulus(4), Delay(2), Response(3)
if ~exist('cInd','var') || isempty(cInd)
    cInd = true(1,length(bhv.Rewarded));
end
if ~exist('dateRange','var') || isempty(dateRange)
    dateRange = [1 inf];
end

cInd = cInd & bhv.date >= dateRange(1) & bhv.date <= dateRange(2);
oInd = cInd & bhv.optoDur == 0.3 & bhv.optoType > 0 & bhv.optoSide == 3 & bhv.DistStim == 0; %bilateral stimulation in detection trials
prevInd = bhv.SessionNr(find(cInd) - 1) == bhv.SessionNr(cInd); %only use trials where the previous trial was frome the same session

for iOutcome = 1:2
    
    bhv.Rewarded
    
    % non-opto performance
    pInd = cInd & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0 & bhv.Rewarded == iOutcome; %only use active detection trials on current side
    dInd = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(oInd))); %non-optogenetic trials from same sessions
    dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
    rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
    
    allData.detect(iOutcome) = sum(rInd)/sum(dInd); %percent right choices
    [allData.detectUp(iOutcome), allData.detectLow(iOutcome)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
    allData.stimCnt(iOutcome) = sum(dInd);
    
    for iTime = 1 : 5 %different task episodes.
        for stimLoc = 1 : 3 % opto performance - frontal (1) / parietal (2) / S1 (3)
            
            dInd = pInd & oInd & bhv.stimLocation == stimLoc & bhv.optoType == typeOrder(iTime); %optogenetic trials
            dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
            rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
            allData.optoDetect(iOutcome,iTime,stimLoc) = sum(rInd)/sum(dInd); %percent right choices
            [allData.optoDetectUp(iOutcome,iTime,stimLoc), allData.optoDetectLow(iOutcome,iTime,stimLoc)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
            
        end
    end
end

% same thing for all trials non-opto performance
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0; %use active detection trials on both sides
dInd = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(oInd))); %non-optogenetic trials from same sessions
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
rInd = dInd & bhv.Rewarded; %correct trials

allData.detect(iOutcome+1) = sum(rInd)/sum(dInd); %percent correct choices
[allData.detectUp(iOutcome+1), allData.detectLow(3)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
allData.stimCnt(iOutcome+1) = sum(dInd);
allData.controlCnt = sum(dInd); 

for iTime = 1 : 5 %different task episodes.
    for stimLoc = 1 : 3 % opto performance - frontal / parietal / S1
        
        dInd = pInd & oInd & bhv.stimLocation == stimLoc & bhv.optoType == typeOrder(iTime); %optogenetic trials
        dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
        rInd = dInd & bhv.Rewarded; %correct trials
        allData.optoDetect(iOutcome+1,iTime,stimLoc) = sum(rInd)/sum(dInd); %percent right choices
        [allData.optoDetectUp(iOutcome+1,iTime,stimLoc), allData.optoDetectLow(iOutcome+1,iTime,stimLoc)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
        allData.optoCnt(iTime,stimLoc) = sum(dInd);
        
    end
end