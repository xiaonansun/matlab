function waveform = buildAudioWaveForm(eventList, eventDuration, samplingFreq)
% white noise sound event
event_part = (rand(1,samplingFreq*eventDuration) * 2) - 1;
% event_part = 2*GenerateSineWave(samplingFreq, 1000, eventDuration);
silence_part = zeros(1,length(event_part));
waveform = [];
for ev = 1:length(eventList)
    this_event = eventList(ev);
    if this_event == 1
        waveform = [waveform event_part silence_part];
    else
        waveform = [waveform silence_part silence_part];
    end
end
end