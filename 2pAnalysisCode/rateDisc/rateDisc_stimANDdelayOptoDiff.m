function [out, earlyOff] = rateDisc_stimANDdelayOptoDiff(bhv)
%short function to compute detection performance with optogenetic
%inactiavtion during stim+delay period and compare performance against
%non-inactivation trials. Returns performance for each session and all
%trials combined.
% bhv is the behavioral array from bpod, cInd indicates trials be used in
% the all-trials-combined analysis.

%% for each animal, determine the early and late stage in training and get behavior for early or late inactivation
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0; %use active detection trials on both sides
oInd = (bhv.optoDur == 1.5 | bhv.optoDur == 1.7) & bhv.optoType == 1 & bhv.optoSide == 3 & bhv.DistStim == 0 & bhv.optoPower > 1; %bilateral stim+delay stimulation.

out.detect = NaN(2,2,3,length(unique(bhv.AnimalID))); % early/late x ctrl/opto x location x animals
out.error = NaN(2,2,2,3,length(unique(bhv.AnimalID))); % early/late x ctrl/opto x up/down x location x animals
for iAnimals = unique(bhv.AnimalID)
    aInd = pInd & ismember(bhv.AnimalID,iAnimals); %performed trials for current animal
    cDates = bhv.date(oInd & aInd);
    
    b = find(diff(cDates) < -10, 2, 'last');  %find time differences that are larger than 20 days
    if ~isempty(b)
        earlyOff = cDates(b(1)+1); %recordings before jump are early
        earlyInd = bhv.date <= earlyOff;
        
        lateOn = cDates(b(end)); %late recordings. sometimes there are two jumps here.
        lateInd = bhv.date >= lateOn; %divide session early and late in training
    elseif isempty(b) && ~isempty(cDates) %no jumps = no late sessions
        earlyInd = bhv.date <= cDates(1);
        lateInd = false(1,length(oInd));
    else
        earlyInd = [];
    end
    
    if ~isempty(earlyInd)
        % get performance for early and late inacivation at different locations. 1 = frontal, 2 = parietal, 3 = S1/V1
        for stimLoc = 1 : 3
            dInd = aInd & bhv.stimLocation == stimLoc; %trials for current location and animal
            rInd = dInd & bhv.Rewarded; %correct trials at current location and animal
            [out.detect(:,:,stimLoc,iAnimals), out.error(:,:,:,stimLoc,iAnimals)] = ...
                rateDisc_earlyAndLate(rInd, dInd, oInd, earlyInd, lateInd);
        end
        
        % same thing for individual sessions
        cSessions = unique(bhv.SessionNr(oInd & aInd)); %opto sessions
        out.sDetect{iAnimals} = NaN(2,length(cSessions)); % ctrl/opto x session
        out.sError{iAnimals} = NaN(2,2,length(cSessions)); % ctl/opto x up/down x session
        out.sInfo{iAnimals} = NaN(length(cSessions),3); %this is the identifier for stimulation location, early or late and trialcount
        Cnt = 0;
        for iSessions = cSessions
            Cnt = Cnt +1;
            
            dInd = pInd & bhv.SessionNr == iSessions; %trials for current session
            rInd = dInd & bhv.Rewarded; %correct trials at current location and animal
            
            out.sInfo{iAnimals}(Cnt,1) = unique(bhv.stimLocation(dInd)); %get current location
            out.sInfo{iAnimals}(Cnt,2) = any(earlyInd(dInd)) + any(lateInd(dInd))*2; %get time of training. 1 = early, 2 = late
            out.sInfo{iAnimals}(Cnt,3) = sum(dInd); %trialcount for current session
            
            [a,b] = rateDisc_earlyAndLate(rInd, dInd, oInd, earlyInd, lateInd); % same analysis as above but only early or late will make sense here
            out.sDetect{iAnimals}(:,Cnt) = nanmean(a,1);
            out.sError{iAnimals}(:,:,Cnt) = squeeze(nanmean(b,1));
        end
    end
end

