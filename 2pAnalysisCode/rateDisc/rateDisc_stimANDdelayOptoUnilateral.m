function allData = rateDisc_stimANDdelayOptoUnilateral(bhv,cInd,dateRange)
%short function to compute detection performance with optogenetic
%inactiavtion during stim+delay period in two cortical locations.

if ~exist('cInd','var') || isempty(cInd)
    cInd = true(1,length(bhv.Rewarded));
end
if ~exist('dateRange','var') || isempty(dateRange)
    dateRange = [1 inf];
end

cInd = cInd & bhv.date >= dateRange(1) & bhv.date <= dateRange(2);
oInd = cInd & bhv.optoDur == 1.5 & bhv.optoType == 1 & bhv.optoSide ~= 3 & bhv.DistStim == 0; %unilateral stim+delay stimuation.

for iSide = 1:2
    % non-opto performance
    pInd = cInd & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0 & bhv.CorrectSide == iSide; %only use active detection trials on current side
    dInd = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(oInd))); %non-optogenetic trials from same sessions
    dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
    rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
    
    allData.detect(iSide) = sum(rInd)/sum(dInd); %percent right choices
    [allData.detectUp(iSide), allData.detectLow(iSide)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
    allData.stimCnt(iSide) = sum(dInd);
    
    for optoSide = 1:3 % optogenetic stimulation side - left / right
        for stimLoc = 1 : 3 % opto performance - frontal / parietal / S1
            if optoSide == 3
                dInd = pInd & oInd & bhv.stimLocation == stimLoc ; %both sides
            else
                dInd = pInd & oInd & bhv.stimLocation == stimLoc & bhv.optoSide == optoSide ; %optogenetic trials
            end
            dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
            rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
            allData.optoDetect(iSide,stimLoc,optoSide) = sum(rInd)/sum(dInd); %percent right choices
            [allData.optoDetectUp(iSide,stimLoc,optoSide), allData.optoDetectLow(iSide,stimLoc,optoSide)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
            
        end
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

for optoSide = 1 : 3 % optogenetic stimulation side - left / right or both
    for stimLoc = 1 : 3 % opto performance - frontal / parietal / S1
        if optoSide == 3
            dInd = pInd & oInd & bhv.stimLocation == stimLoc; %optogenetic trials
        else
            dInd = pInd & oInd & bhv.stimLocation == stimLoc & bhv.optoSide == optoSide; %optogenetic trials
        end
        dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
        rInd = dInd & bhv.Rewarded; %correct trials
        allData.optoDetect(iSide+1,stimLoc,optoSide) = sum(rInd)/sum(dInd); %percent correct choices
        [allData.optoDetectUp(iSide+1,stimLoc,optoSide), allData.optoDetectLow(iSide+1,stimLoc,optoSide)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
        allData.optoCnt(stimLoc,optoSide) = sum(dInd);
    end
end

% get optogenetic effects for ipsi- and contraleteral side
for stimLoc = 1 : 3 % opto performance - frontal / parietal / S1
    
    pInd = cInd & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0; %only use active detection trials on current side
    
    % do comparison for all trials or only correct or incorrect trials
    for x = 1 : 3
        
        dInd = pInd & oInd & bhv.stimLocation == stimLoc; %all optogenetic trials
        dInd = rateDisc_equalizeTrialsPerMouse(bhv, dInd); %equalize trial counts for animals in current selection
        if x == 2
            dInd = dInd & bhv.Rewarded; %only rewarded trials
        elseif x == 3
            dInd = dInd & ~bhv.Rewarded; %only error trials
        end
        
        dInd = rateDisc_equalizeTrials(dInd, bhv.CorrectSide == 1, bhv.optoSide == 1); %equalize left and right target stimuli

        % amount of choices, contralateral to stimulation side
        rInd = dInd & bhv.optoSide ~= bhv.ResponseSide; %contralateral choice trials
        allData.optoContra(stimLoc,x) = sum(rInd)/sum(dInd); %percent contra choices
        allData.optoContraCnt(stimLoc,x) = sum(dInd); %number of trials used to get contra choice behavior
        [allData.optoContraUp(stimLoc,x), allData.optoContraLow(stimLoc,x)] = Behavior_wilsonError(sum(rInd), sum(dInd)); %error
    end
end