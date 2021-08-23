
keep SessionData
for iTrials = 1 : length(SessionData.Rewarded)
    
    trialStart(iTrials) = SessionData.RawEvents.Trial{iTrials}.States.PlayStimulus(1);
    trialDur(iTrials) = SessionData.RawEvents.Trial{iTrials}.States.StopStim(end);
    
    varStimOn(iTrials) = SessionData.RawEvents.Trial{iTrials}.States.RunningStimulus(1) - SessionData.RawEvents.Trial{iTrials}.States.PlayStimulus(1);
    stimDur(iTrials) = diff(SessionData.RawEvents.Trial{iTrials}.States.RunningStimulus);
    delayDur(iTrials) = diff(SessionData.RawEvents.Trial{iTrials}.States.DecisionWait);

    stimTrialDur(iTrials) = SessionData.RawEvents.Trial{iTrials}.States.PrepareNextTrial(1) - SessionData.RawEvents.Trial{iTrials}.States.RunningStimulus(1);

end
disp(sum(~SessionData.DidNotLever & ~SessionData.DidNotChoose))