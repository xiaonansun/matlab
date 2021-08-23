function bhv = rateDisc_makeBlockBhv(cPath, animal, trainingRange)

if isempty('trainingRange') || isempty(trainingRange)
    trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
end

if ~strcmpi(cPath(end), filesep)
    cPath = [cPath filesep];
end

fPath = [cPath animal filesep 'blockData' filesep];
load([fPath 'trialInfo_' trainingRange '.mat'],'bhvTrials','recs')

bhv = [];
for iRecs = 1 : length(recs)
    
    fPath = [cPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep];
    bhvFile = dir([fPath filesep animal '*.mat']);
    
    load([bhvFile.folder filesep bhvFile.name])
    SessionData.Notes = SessionData.Notes(1:length(SessionData.Rewarded));
    SessionData.MarkerCodes = SessionData.MarkerCodes(1:length(SessionData.Rewarded));
    
    SessionData = selectBehaviorTrials(SessionData,bhvTrials{iRecs});
    SessionData.SessionNr = repmat(iRecs,size(SessionData.Rewarded));
    bhv = appendBehavior(bhv,SessionData); %append into larger array

end

% save back to blockdata path
fPath = [cPath animal filesep 'blockData' filesep];
save([fPath 'blockBhv_' trainingRange '.mat'],'bhv');
