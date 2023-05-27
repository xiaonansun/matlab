function E = behavior_getStimEvents(bhv)
%% Input: bhv is the Bpod behavior session structure (e.g. SessionData)

events = bhv.stimEvents;

%%
% Auditory event timing. Each cell is a vector of event timing relative to
% the onset of the stimulus epoch
% column 1: left-sided event timing while column 2 is right-sided event
% timing
% rows are trials
eventTiming = cellfun(@(x) x{:},[cellfun(@(x) x(1),events,'UniformOutput',false)' cellfun(@(x) x(2),events,'UniformOutput',false)'],'UniformOutput',false); 

% The of count the number of auditory events emerging from the left or right speakers and report them in two columns (column 1 for left and 2 for right)
numEvents= cellfun(@numel,eventTiming); 

% ratio of left to all events presented during a trial; 0=100% LEFT-SIDED events, 1=100% RIGHT-SIDED events
ratioEvents= numEvents(:,2)./sum(numEvents,2); 

E.eventTiming = eventTiming;
E.numEvents = numEvents;
E.ratioEvents = ratioEvents;