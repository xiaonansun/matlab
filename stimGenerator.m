clear all
n = 10000;

current_rate = 5;
pulse_width = 20;
stim_duration = 1000;
interpulse_int = 0;

for i = 1:n
    [waveform_bin,peak_idx] = base_waveform_generator(current_rate,pulse_width,stim_duration);
    num_of_peaks(i) = numel(peak_idx);
    interpulse_int = [interpulse_int diff(peak_idx)];
end
interpulse_int = interpulse_int(2:end);

subplot(2,1,1)
histogram(num_of_peaks,100,'Normalization','probability');
hold on;
subplot(2,1,2)
histogram(interpulse_int,100,'Normalization','probability')