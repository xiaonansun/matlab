function [Performance,bhv] = rateDisc_sessionPerformance(animal,cPath,recs)
% Analze behavioral data from rate discrimination paradigm to get basic readout of animals performance.
% Gets behavioral data from recordings in cPath. This is to get behavioral performance during widefield imaging.

modId = [2 4 6]; %id for modalities. Default is audio, somatosensory and audiosomatosensory
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings; %basic training info
Performance.recType = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs);
bhv = [];

% get data from individual recordings
for iRecs = 1 : length(recs)
    
    fPath = [cPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep];
    bhvFile = dir([fPath animal '*SpatialDisc*.mat']);
    load([fPath bhvFile(1).name], 'SessionData'); %load current bhv file

    Performance.AllTrials(iRecs) = length(SessionData.Rewarded);
    Performance.doneTrials(iRecs) = sum(~SessionData.DidNotChoose & ~SessionData.DidNotLever & logical(SessionData.Assisted));
    Performance.date(iRecs) = datenum(bhvFile(1).name(length([animal '_SpatialDisc_'])+1:length([animal '_SpatialDisc_'])+10)); %isolate dates from Filenames

    for iMod = 1 : 3
        %% get some single session performance data
        ind = ~SessionData.DidNotChoose & ~SessionData.DidNotLever & logical(SessionData.Assisted) & SessionData.StimType == modId(iMod); %only use active trials
        Performance.SelfPerformed(iMod,iRecs) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        
        dInd = logical(SessionData.Assisted & SessionData.DistStim == 0); %index for trials that were detection only
        Performance.Detection(iMod,iRecs) = sum(SessionData.Rewarded(ind & dInd))/sum(SessionData.Rewarded(ind & dInd)+SessionData.Punished(ind & dInd)); %peformance for detection trials
        Performance.Discrimination(iMod,iRecs) = sum(SessionData.Rewarded(ind & ~dInd))/sum(SessionData.Rewarded(ind & ~dInd)+SessionData.Punished(ind & ~dInd)); %peformance for discrimination trials
        Performance.DetTrials(iMod,iRecs) = sum(ind & dInd); %number of detection trials
        Performance.DiscTrials(iMod,iRecs) = sum(ind & ~dInd); %number of discrimination trials
        
        lInd = SessionData.CorrectSide == 1; %index for left-choice trials
        Performance.LeftPerformed(iMod,iRecs) = sum(SessionData.Rewarded(ind & lInd))/sum(SessionData.Rewarded(ind & lInd)+SessionData.Punished(ind & lInd));
        Performance.RightPerformed(iMod,iRecs) = sum(SessionData.Rewarded(ind & ~lInd))/sum(SessionData.Rewarded(ind & ~lInd)+SessionData.Punished(ind & ~lInd));
    end
    
    %% combine into one larger array
    SessionData.SessionNr = repmat(iRecs,1,SessionData.nTrials); %tag all trials in current dataset with session nr
    bhv = appendBehavior(bhv,SessionData); %append into larger array
        
end