% PsychToolboxSoundServer('init'); % initialize once
% cannot use channels 3-6 unless PsychToolboxSoundServer is modified
% Channel 1: audio
% Channel 2: visual

thisRate = 5;
show_audio = 0;
show_visual = 1;
visual_rate = 20;
audio_rate = 5;
eventDuration = 0.020;

difficulty = 1;
min_rate = 5;
max_rate = 20;
stim_duration= 1;
event_bin_duration = 0.040;
event_duration = 0.020;
sampling_freq = 192000;

pwm = 1;
brightness = 20;
soundLoudness = 1;


% Build Visual Waveforms
% [eventList_viusual_right,got_it]= getRandomStimEvents(visual_rate_right, stimDuration, eventDuration);
% [eventList_viusual_left,got_it]= getRandomStimEvents(visual_rate_left, stimDuration, eventDuration);
% visual_waveform_right = buildVisualWaveForm(eventList_viusual_right, eventDuration, samplingFreq, brightness, pwm);
% visual_waveform_left = buildVisualWaveForm(eventList_viusual_left, eventDuration, samplingFreq, brightness, pwm);
% visual_waveform = [visual_waveform_right;visual_waveform_left];
% [eventList_viusual,got_it]= getRandomStimEvents(visual_rate, stimDuration, eventDuration);
% visual_waveform = buildVisualWaveForm(eventList_viusual, eventDuration, samplingFreq, brightness, pwm);
stim_events = generate_stim_events(difficulty, event_bin_duration, stim_duration, min_rate, max_rate);
visual_waveform = create_visual_waveform(stim_events, sampling_freq, brightness, event_duration, stim_duration);

% Build Audio Waveforms
% [eventList_audio_right,got_it]= getRandomStimEvents(visual_rate_right, stimDuration, eventDuration);
% [eventList_audio_left,got_it]= getRandomStimEvents(visual_rate_left, stimDuration, eventDuration);
% audio_waveform_right = buildAudioWaveForm(eventList_audio_right, eventDuration, samplingFreq);
% audio_waveform_left = buildAudioWaveForm(eventList_audio_left, eventDuration, samplingFreq);
% audio_waveform = [audio_waveform_right;audio_waveform_left];

[eventList_audio,got_it]= getRandomStimEvents(audio_rate, stimDuration, eventDuration);
audio_waveform = buildAudioWaveForm(eventList_audio, eventDuration, samplingFreq);

if length(audio_waveform) < length(visual_waveform)
    visual_waveform = visual_waveform(:,1:length(audio_waveform));
end
waveform = [show_visual*visual_waveform;show_audio*audio_waveform];

time = 1/samplingFreq:1/samplingFreq:length(waveform)/samplingFreq;

plot(time,waveform(1,:),'-b'); hold on;
plot(time,waveform(2,:),'-r');
%% upload to PsychToolboxSoundServer

PsychToolboxSoundServer('load',1,waveform);

%% play sound loaded in previous step
PsychToolboxSoundServer('play',1);

%% Test start trial sound
% 14100 Hz is audible limit

% waveStartSound = 0.1 * soundLoudness * GenerateSineWave(samplingFreq, 7000, 0.1);
% 
% PsychToolboxSoundServer('load',1,waveStartSound);
% PsychToolboxSoundServer('play',1);
% 
