function allData = rateDisc_stimANDdelayOpto(bhv,cInd,dateRange)
%short function to compute detection performance with optogenetic
%inactiavtion during stim+delay period in two cortical locations.

if ~exist('cInd','var') || isempty(cInd)
    cInd = true(1,length(bhv.Rewarded));
end
if ~exist('dateRange','var') || isempty(dateRange)
    dateRange = [1 inf];
end

cInd = cInd & bhv.date >= dateRange(1) & bhv.date <= dateRange(2);
oInd = cInd & bhv.optoDur == 1.5 & bhv.optoType == 1 & bhv.optoSide == 3 & bhv.DistStim == 0; %bilateral stim+delay stimuation.
for iSide = 1:2
    % non-opto performance
    pInd = cInd & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0 & bhv.CorrectSide == iSide; %only use active detection trials on current side
    dInd = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(oInd))); %non-optogenetic trials from same sessions
    dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
    rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
    
    allData.detect(iSide) = sum(rInd)/sum(dInd); %percent right choices
    [allData.detectUp(iSide), allData.detectLow(iSide)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
    allData.stimCnt(iSide) = sum(dInd);
    
    % opto performance - frontal / parietal / S1
    for stimLoc = 1 : 3
        dInd = pInd & oInd & bhv.stimLocation == stimLoc; %optogenetic trials
        dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
        rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
        allData.optoDetect(iSide,stimLoc) = sum(rInd)/sum(dInd); %percent right choices
        [allData.optoDetectUp(iSide,stimLoc), allData.optoDetectLow(iSide,stimLoc)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
    end
end

% same thing for all trials non-opto performance
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0; %use active detection trials on both sides
dInd = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(oInd))); %non-optogenetic trials from same sessions
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
rInd = dInd & bhv.Rewarded; %correct trials

allData.detect(iSide+1) = sum(rInd)/sum(dInd); %percent correct choices
[allData.detectUp(iSide+1), allData.detectLow(3)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
allData.stimCnt(iSide+1) = sum(dInd);
allData.optoCnt(stimLoc +1) = sum(dInd); 

% opto performance - frontal / parietal / S1
for stimLoc = 1 : 3
    dInd = pInd & oInd & bhv.stimLocation == stimLoc; %optogenetic trials
    dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
    rInd = dInd & bhv.Rewarded; %correct trials
    allData.optoDetect(iSide+1,stimLoc) = sum(rInd)/sum(dInd); %percent right choices
    [allData.optoDetectUp(iSide+1,stimLoc), allData.optoDetectLow(iSide+1,stimLoc)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
    allData.optoCnt(stimLoc) = sum(dInd);
end