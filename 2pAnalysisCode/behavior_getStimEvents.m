function E = behavior_getStimEvents(bhv)
%% Input: bhv is the Bpod behavior session structure (e.g. SessionData)


%%
% Auditory event timing. Each cell is a vector of event timing relative to
% the onset of the stimulus epoch
% Column 1: left-sided event timing 
% column 2: right-sided event timing
% Each row is a single trials
events = bhv.stimEvents;
eventTiming = cellfun(@(x) x{:},[cellfun(@(x) x(1),events,'UniformOutput',false)' cellfun(@(x) x(2),events,'UniformOutput',false)'],'UniformOutput',false); 

% The of count the number of auditory events emerging from the left or right speakers and report them in two columns (column 1 for left and 2 for right)
numEvents= cellfun(@numel,eventTiming); 

% Proportion of left-sided events
% 0: All events are RIGHT-SIDED
% 1: All events are LEFT-SIDED 
propLeftEvents= numEvents(:,1)./sum(numEvents,2); 

E.eventTiming = eventTiming;
E.numEvents = numEvents;
E.ratioEvents = propLeftEvents;