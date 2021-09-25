function rTimes = Behavior_reactionTimes(SessionData)
% compute reaction times

rTimes = zeros(1,SessionData.nTrials);

for iTrials = 1:SessionData.nTrials
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lLick = SessionData.RawEvents.Trial{iTrials}.Events.Port1In; %left licks
    else
        lLick = []; %no licks
    end
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port3In') %check for licks
        rLick = SessionData.RawEvents.Trial{iTrials}.Events.Port3In; %right licks
    else
        rLick = []; %no licks
    end
    
    if (isempty(lLick) && isempty(rLick)) || any(any(isnan(SessionData.RawEvents.Trial{iTrials}.States.WaitForResponse)))
        rTimes(:,iTrials) = NaN; %no response
    else
        window = SessionData.RawEvents.Trial{iTrials}.States.WaitForResponse;
        lLick = lLick-window(1); lLick(lLick < 0 | lLick > (window(2) - window(1))) = []; %compute left response during decision window
        if isempty(lLick); lLick = inf; end
        rLick = rLick-window(1); rLick(rLick < 0 | rLick > (window(2) - window(1))) = []; %compute right response during decision window
        if isempty(rLick); rLick = inf; end

        if min(lLick) <  min(rLick)  % if left licks are first
            rTimes(iTrials) = lLick(1);
        elseif min(lLick) >  min(rLick)   % if right licks are first
            rTimes(iTrials) = rLick(1);
        else
            rTimes(iTrials) = NaN; %no response
        end
        clear window
    end
    clear  lLick rLick
end