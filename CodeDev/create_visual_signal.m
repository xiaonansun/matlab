function signal_visual = create_visual_signal(isis, rate, brightness, flash_duration, long_isi, short_isi)

timevec = (0 : 1/rate : flash_duration);
freq = 200*brightness; % we can't detect a 300 Hz flicker
sine_wave = sin(2*pi * freq * timevec);

[pks,locs] = findpeaks(sine_wave);

flash = zeros(size(sine_wave));
flash(locs) = 1;

%     shift = brightness .* 2 - 1;
%     timevec = (0 : 1/rate : flash_duration);
%     freq = 10000; % we can't detect a 300 Hz flicker
%     sine_wave = sin(2*pi * freq * timevec);
%
%     on_off = sine_wave > - shift;
%     %flash = sine_wave .* on_off;
%     flash = on_off;
%
%     %flash = sine_wave;
%     % figure; plot(flash);

long_dark = zeros(1, round(rate * long_isi));
short_dark = zeros(1, round(rate * short_isi));

signal_visual = flash;
for i = isis
    if i == 1
        signal_visual = [signal_visual short_dark flash];
    elseif i == 2
        signal_visual = [signal_visual long_dark flash];
    end
end

end