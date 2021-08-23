function allData = rateDisc_stimORdelayOpto(bhv,cInd)
%short function to compute detection performance with optogenetic
%inactiavtion during stim+delay period in two cortical locations.

if isempty(cInd)
    cInd = ones(1,length(bhv.Rewarded));
end

oInd = cInd & round(bhv.optoDur,4) >= 0.6 & round(bhv.optoDur,4) <= 0.8 & bhv.optoSide == 3 & bhv.DistStim == 0; %bilateral stim or delay stimuation.
pInd = cInd & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0; %only use active detection trials on current side
dInd = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(oInd))); %non-optogenetic trials from same sessions
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(1) = sum(dInd);
allData.detect = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.detectUp, allData.detectLow] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error

dInd = pInd & oInd & bhv.stimLocation == 1 & bhv.optoType == 1; %optogenetic; frontal stim only
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(2) = sum(dInd);
allData.opto(1,1) = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.optoUp(1,1), allData.optoLow(1,1)] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error

dInd = pInd & oInd & bhv.stimLocation == 1 & bhv.optoType == 2; %optogenetic; frontal delay only
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(3) = sum(dInd);
allData.opto(1,2) = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.optoUp(1,2), allData.optoLow(1,2)] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error

dInd = pInd & oInd & bhv.stimLocation == 2 & bhv.optoType == 1; %optogenetic; parietal stim only
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(4) = sum(dInd);
allData.opto(2,1) = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.optoUp(2,1), allData.optoLow(2,1)] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error

dInd = pInd & oInd & bhv.stimLocation == 2 & bhv.optoType == 2; %optogenetic; parietal delay only
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(5) = sum(dInd);
allData.opto(2,2) = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.optoUp(2,2), allData.optoLow(2,2)] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error

dInd = pInd & oInd & bhv.stimLocation == 3 & bhv.optoType == 1; %optogenetic; S1 stim only
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(6) = sum(dInd);
allData.opto(3,1) = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.optoUp(3,1), allData.optoLow(3,1)] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error

dInd = pInd & oInd & bhv.stimLocation == 3 & bhv.optoType == 2; %optogenetic; S1 delay only
dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
allData.tCounts(7) = sum(dInd);
allData.opto(3,2) = sum(bhv.Rewarded(dInd))/sum(dInd); %percent correct choices
[allData.optoUp(3,2), allData.optoLow(3,2)] = Behavior_wilsonError(sum(bhv.Rewarded(dInd)), sum(dInd)); %error