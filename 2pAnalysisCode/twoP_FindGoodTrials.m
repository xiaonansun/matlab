function [trials, nTrials] = twoP_FindGoodTrials(cPath)
% Code to identify trials that should be used for imaging analysis. 
% Reject incomplete trials.
       
%% load behavior data, check stim delay and add indicies for potentially unfit baseline activity
bhvFile = strsplit(cPath,filesep);
bhvFile = dir([cPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
load([cPath bhvFile.name]); %load behavior data

badTrial = false(1,SessionData.nTrials);
noStimTrial = false(1,SessionData.nTrials);
SessionData.TouchCnt = zeros(1,length(SessionData.Rewarded));
trials = 1 : length(SessionData.Rewarded);
for iTrials = 1 : length(trials)
    
    stimDelay = [];
    stimTime = SessionData.RawEvents.Trial{iTrials}.States.RunningStimulus;
    % check for delayed stimulus responses
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Wire3High')
        stimDelay = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High - stimTime(1); %recorded stimOn minus presumed stimOn
    end
    
    stimDelay(stimDelay < 0) = [];
    if any(stimDelay > 0.1) %delay above 100ms is not ok
        badTrial(iTrials) = true;
    end
    
    if  isempty(stimDelay) % no stimulus
        noStimTrial(iTrials) = true;
    end
    
    leverTouch = [];
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Wire1Low')
        leverTouch = SessionData.RawEvents.Trial{iTrials}.Events.Wire1Low;
    end
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Wire2Low')
        leverTouch = [leverTouch SessionData.RawEvents.Trial{iTrials}.Events.Wire2Low];
    end
    SessionData.TouchCnt(iTrials) = sum(((sort(leverTouch) - stimTime(1))) < 0 &  ((sort(leverTouch) - stimTime(1))) >= -(SessionData.Settings.preStimDelay + 1));
end

if any(badTrial)
    fprintf('Found %d/%d bad trials due to stimOn delays above 100ms.\n',sum(badTrial),length(badTrial))
    trials(badTrial) = [];
end

if any(noStimTrial)
    fprintf('Removed %d/%d bad trials due to missing stimulus onset.\n',sum(noStimTrial),length(noStimTrial))
    trials(noStimTrial) = [];
end

save([cPath bhvFile.name], 'SessionData'); %save behavior data
nTrials = SessionData.nTrials;