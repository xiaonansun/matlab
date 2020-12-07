function waveform = buildVisualWaveForm(eventList, eventDuration, samplingFreq, brightness, pwm)
%07-Sept-2016: added pwm flag. for new ALA Scientific LED panel driver

% visual event
event_part = GenerateSineWave(samplingFreq, 200*brightness, eventDuration);

if pwm
    [pks,locs] = findpeaks(event_part);
    
    flash = zeros(size(event_part));
    flash(locs) = 1;
    event_part = flash;
else
    event_part(event_part < 0.975) = 0;
end

if numel(eventList) == 1
    waveform = event_part;
    return
end

silence_part = zeros(1,length(event_part));
waveform = [];
for ev = 1:length(eventList)
    this_event = eventList(ev);
    if this_event == 1
        waveform = cat(2,waveform, event_part, silence_part);
    else
        waveform = cat(2,waveform, silence_part, silence_part);
    end
end
end