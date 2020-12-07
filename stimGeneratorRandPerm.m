clear all;
stim_duration = 1000;
num_of_pulses = 5;
pulse_duration = 20;
interpulse_int = 20;
min_ipi = pulse_duration+interpulse_int;

ipi_all_trials = 0;
for i = 1:10000
ipi=0;
while sum(ipi < min_ipi) >= 1
stim = [ones(1,num_of_pulses-1) zeros(1,stim_duration-num_of_pulses)];
stim = [1 stim(randperm(length(stim)))];
ipi=diff(find(stim==1));
end
ipi_all_trials=[ipi_all_trials ipi];
ipi_all_trials=ipi_all_trials(2:end);
end

%% plot distribution
[muhat,muci]=expfit(ipi_all_trials);
h1=histogram(ipi_all_trials,1000,'Normalization','pdf');
h1.FaceColor='none';

fig1=gcf;
ax1=fig1.CurrentAxes;
ax1.Title.String=['Interpulse interval for ' num2str(num_of_pulses) ' pulses'];
ax1.XLabel.String='Time (ms)'; ax1.YLabel.String='Probability density';
ax1.FontSize=12; ax1.FontName='Arial';
ax1.TickDir='out'; ax1.Box='off'; ax1.LineWidth=1; ax1.Clipping='off'; ax1.Color='none';
ax1.OuterPosition=[0 0 1 1];
