function bhv = rateDisc_collectBehaviorSession(dPath, recs, animal, bhvTrials)

bhv = [];
for iRecs = 1 : length(recs)
    
    fPath = [dPath recs(iRecs).name filesep];
    bhvFile = dir([fPath filesep animal '*.mat']);
    
    load([bhvFile.folder filesep bhvFile.name])
    SessionData.Notes = SessionData.Notes(1:length(SessionData.Rewarded));
    SessionData.MarkerCodes = SessionData.MarkerCodes(1:length(SessionData.Rewarded));
    
    SessionData = selectBehaviorTrials(SessionData,bhvTrials{iRecs});
    SessionData.SessionNr = repmat(iRecs,size(SessionData.Rewarded));
    bhv = appendBehavior(bhv,SessionData); %append into larger array

end