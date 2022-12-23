function [trigger_idx,aligned_sync_events,sync_new] = np_extractSyncEvents(sync, sampling_rate)

%%
thresh = 5000;
min_dur = 4; % minimum duration in seconds between each trigger event
interval = sampling_rate*min_dur;

idxDelta = [0;find(abs(diff(sync)) > thresh);length(sync)];
% z = sync(1:idxDelta(1))- max(sync(1:idxDelta(1)));
% y = sync(idxDelta(1)+1:idxDelta(2))- max(sync(idxDelta(1)+1:idxDelta(2)));

% The following is one way of dealing with baseline changing multiple
% times.
% (1) Find the large, instantaneous changes in signals as specified by
% 'thresh'
% (2) Take each consecutive intervals bounded by these instantaneous
% changes and subtract that by the maximum value of that interval (hence normalizing each interval to a max of 0)

z = [];
for i = 1:length(idxDelta)-1
    z = [z;sync(idxDelta(i)+1:idxDelta(i+1)) - max(sync(idxDelta(i)+1:idxDelta(i+1)))]; % subtracts every interval by the max value of that interval
end

x = ~ismember(z,unique(z(1:30000))); %The oscillatory signals only occur in steps of 4. View this with plot(z(1:30000)).
% So taking the first second of oscillatory values wouldn't capture any event signals. 
% Therefore, if any signal doesn't match these unique oscillatory values is considered the stimulus signal.
sync_new = x;

sync_new = x;
trigger_idx = [];
cnt = 1;
while cnt < length(sync_new)
    if sync_new(cnt) == 1
        trigger_idx = [trigger_idx cnt];
        cnt = cnt+interval;
    else
        cnt = cnt+1;
    end
end

aligned_sync_window = [0 0.05]; % window for aligning sync events into a matrix
aligned_sync_events = zeros(length(trigger_idx),abs(aligned_sync_window(1)*sampling_rate-1)+aligned_sync_window(2)*sampling_rate);

for i = 1:length(trigger_idx)
   aligned_sync_events(i,:) = 1-[sync_new(trigger_idx(i)+aligned_sync_window(1)*sampling_rate:trigger_idx(i)-1,1)' sync_new(trigger_idx(i):trigger_idx(i)+aligned_sync_window(2)*sampling_rate,1)'];
end

% syncA = [0;sync];
% dSyncA = diff(syncA);
% idxDeltaA = [0;find(abs(dSyncA) > thresh);length(syncA)];
% 
% for i = length()