function [eventList,got_it]= getRandomStimEvents(thisRate, stimDuration, eventDuration)
frames = floor(stimDuration/eventDuration/2);
maxAttempts = 1000;
got_it = 0;
for attempt = 1:maxAttempts
    eventList = [1 (poissrnd(rand, 1, frames-1) > 0)];
    if sum(eventList) == thisRate
        got_it = 1;
        break;
    end
end

if ~got_it
    %try again to generate the stimulus, using higher lambda value.
    got_it = 0;
    for attempt = 1:maxAttempts
        eventList = [1 (poissrnd(5*rand, 1, frames-1) > 0)];
        if sum(eventList) == thisRate
            got_it = 1;
            break;
        end
    end
end