function [waveform_bin,peak_idx] = base_waveform_generator(current_rate,pulse_width,stim_duration)
% current_rate = 5;
% pulse_width = 20;
% stim_duration = 1000;
peak_idx = 1;
stim = [ones(1,pulse_width) zeros(1,stim_duration-pulse_width)];

j = 2;
for i = 1:stim_duration
    if stim(i) == 0 && sum(stim(i-pulse_width:i)) == 0 && i <= stim_duration-pulse_width
        stim(i) = rand<current_rate/stim_duration;
        if stim(i) == 1
            stim(i:i+pulse_width-1) = ones(1,pulse_width);
            peak_idx(j) = i;
            j=j+1;
        end
    end
end

waveform_bin = stim;
% plot(stim)
% peak_idx