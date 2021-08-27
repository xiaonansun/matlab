function E = behavior_getStimEvents(bhv)

events =bhv.stimEvents;
% organize stimulus events of each trial (row) into columns: column 1 is
% left and column 2 is right
events = cellfun(@(x) x{:},[cellfun(@(x) x(1),events,'UniformOutput',false)' cellfun(@(x) x(2),events,'UniformOutput',false)'],'UniformOutput',false); 
% Auditory event timing. Each cell is a vector of event timing relative to
% the onset of the stimulus epoch

numEvents= cellfun(@numel,events); 
% The of count the number of auditory events emerging from the left or right speakers and report them in two columns (column 1 for left and 2 for right)

ratioEvents= numEvents(:,2)./sum(numEvents,2); 
% ratio of left to all events presented during a trial; 0=100% LEFT-SIDED events, 1=100% RIGHT-SIDED events

E.events = events;
E.numEvents = numEvents;
E.ratioEvents = ratioEvents;