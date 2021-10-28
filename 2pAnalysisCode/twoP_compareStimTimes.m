function twoP_compareStimTimes(D,bhv)
%%
% This function computes the inter-stimulus interval between adjacent
% trials and compares the timing between MScan analog timestamps versus Bpod
% timestamps
unityX = [0 20];
unityY = [0 20];

stimTimeStamp = nan(1,length(bhv.RawEvents.Trial));
for i = 1:length(bhv.RawEvents.Trial)
    try
    stimTimeStamp(i) = bhv.RawEvents.Trial{i}.Events.Wire3High;
    catch
    end
end

trialStimTime = bhv.TrialStartTimestamp+stimTimeStamp;
d_trialStimTime = diff(trialStimTime);
shuf_d_trialStimTime = d_trialStimTime(randperm(length(d_trialStimTime)));
shift_d_trialStimTime = circshift(d_trialStimTime,1);

mscanStimTime = nan(1,length(D.trialNumbers));
mscanStimTime = D.stimSamplesOrig/1000;
d_mscanStimTime = diff(mscanStimTime);

v1 = [unityX(1) unityY(1)];
v2 = [unityX(2) unityY(2)];

for i = 1:length(d_trialStimTime)
ptUnshuffled = [d_trialStimTime(i) d_mscanStimTime(i)];
ptShuffled = [shuf_d_trialStimTime(i) d_mscanStimTime(i)];
ptShifted = [shift_d_trialStimTime(i) d_mscanStimTime(i)];
dUnshuffled(i) = point_to_line_distance(ptUnshuffled,v1,v2);
dShuffled(i) = point_to_line_distance(ptShuffled,v1,v2);
dShifted(i) = point_to_line_distance(ptShifted,v1,v2);
end

%%
fISI = figure(1);
fISI.Position=[500 500 560 420];
plot([unityX(1) unityX(2)],[unityY(1) unityY(2)],...
    'linewidth',4,...
    'Color',[0 0 0]+0.7); hold on; grid on;
plot(shift_d_trialStimTime,d_mscanStimTime,'.b',...
    'MarkerSize',8); 
plot(shuf_d_trialStimTime,d_mscanStimTime,'.r',...
    'MarkerSize',8); 
plot(d_trialStimTime,d_mscanStimTime,'.k',...
    'MarkerSize',8); 
xlabel('Bpod ISI (sec)'); ylabel('2P ISI (sec)');
xlim([0 20]); ylim([0 20]);
legend({'unity','shuffled','shifted','unshuffled'},...
    'Location','Northwest',...
    'Box','off');
ax = gca;
fig_configAxis(ax);
title(['Inter-stimulus interval across consecutive trials: ' D.animal ' ' D.session]);

% Subpanel
axes('Position',[0.85 0.2 .05 .4])
hold on;
marker_size = 2;
scatter(rand(1,length(dShifted)),dShifted,marker_size,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor','b');
scatter(rand(1,length(dShuffled)),dShuffled,marker_size,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor','r');
scatter(rand(1,length(dUnshuffled)),dUnshuffled,marker_size,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor','k'); 
% yline(10e-5);
% ylim([-0.1 20]); 
ax = gca; fig_configAxis(ax); ylabel('Distance from unity');
set(gca, 'YScale', 'log',...
    'Box','off',...
    'XTick',[],...
    'Color','none');
ax.XAxis.Visible = 'off';
% end subpanel

exportgraphics(fISI,fullfile(D.suite2pDir,'inter_stimulus_interval.pdf'));