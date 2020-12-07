%%
plot(sync_b(1:1000000,1))

%%
edges = [1:imec.fsAP*bin_size:imec.nSamplesAP imec.nSamplesAP];
[N,edges]=histcounts(spikeTimes,edges);
trigger_idx_ds = uint32(trigger_idx/(imec.fsAP*bin_size));
trigger_idx_ds = trigger_idx_ds(2:end);

%%
num_of_clusters = length(unique(spike_clusters));

psth_trial = zeros(num_of_clusters,length(time),length(trigger_idx_ds));

for i = 1:num_of_clusters
%    idx = find(spike_clusters == spike_clusters(i));
   spikeTimes_mat{i} = spikeTimes(spike_clusters == spike_clusters(i));
   [N_ds(i,:),edges]=histcounts(spikeTimes_mat{i},edges); % iCluster x activity, iCluster is the index of the cluster ID and activity is the binned neural activity
   
   for j = 1:length(trigger_idx_ds)
       psth_trial(i,:,j) = [N_ds(i,trigger_idx_ds(j)-preStimD/bin_size:trigger_idx_ds(j)-1) N_ds(i,trigger_idx_ds(j):trigger_idx_ds(j)+stimD/bin_size)];
   end
%    display(num2str(i));
end

%% Plotting
time_zero_idx = preStimD/bin_size+1;
xtick_idx = [1 time_zero_idx time_zero_idx+1/bin_size size(psth_trial,2)];
xTickLabel = string(time);
cluster_id = 694;
trial_id = 110;

figure(1);
psth_trial_mean = mean(psth_trial,3); % mean across all trials
int_mean_trial_mean = [mean(psth_trial_mean(:,time_zero_idx:time_zero_idx+1/bin_size),2) (1:1:length(psth_trial_mean))'];
idx_sort_trial_mean = sortrows(int_mean_trial_mean,1);

imagesc(psth_trial_mean(idx_sort_trial_mean(:,2),:)); xlabel('Time (seconds); 0 = Stim onset'); ylabel('Cluster ID');
title('Mean across all trials');
set(gca,'XTick',xtick_idx,...
    'XTickLabel',xTickLabel(xtick_idx));

figure(2);
psth_cluster_mean = mean(psth_trial,1); % mean across all clusters
for i = 1:size(psth_trial,3)
psth_cluster_mean_mat(i,:) = psth_cluster_mean(1,:,i);
end

int_mean = [mean(psth_cluster_mean_mat(:,time_zero_idx:time_zero_idx+1/bin_size),2) (1:1:length(psth_cluster_mean_mat))'];
idx_sort = sortrows(int_mean,1);

imagesc(psth_cluster_mean_mat(idx_sort(:,2),:)); xlabel('Time (seconds); 0 = Stim onset'); ylabel('Trial #');
title('Mean across all clusters');
set(gca,'XTick',xtick_idx,...
    'XTickLabel',xTickLabel(xtick_idx));

figure(3);
imagesc(psth_trial(:,:,trial_id)); xlabel('Time (seconds); 0 = Stim onset'); ylabel('Cluster ID');
title(['Trial ' num2str(trial_id) ', all clusters']);
set(gca,'XTick',xtick_idx,...
    'XTickLabel',xTickLabel(xtick_idx));

figure(4);
for i = 1:size(psth_trial,3)
psth_single_cluster(i,:) = psth_trial(cluster_id,:,i);
end

imagesc(psth_single_cluster); xlabel('Time (seconds); 0 = Stim onset'); ylabel('Trial #');
title(['Cluster ' num2str(cluster_id)]);
set(gca,'XTick',xtick_idx,...
    'XTickLabel',xTickLabel(xtick_idx));

%% Plot PSTHs using spikes

analysis_window = [-1 7]; % look at spike times from 0.3 sec before each event to 1 sec after

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

eventTimes = trigger_idx(2:end)/30000;
depthBinSize = 80; % in units of the channel coordinates, in this case µm
timeBinSize = 0.01; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

% [timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
%     depthBinSize, timeBinSize, eventTimes, analysis_window, bslWin);

trialGroups = ones(size(eventTimes)); 
psthViewer(sp.st, sp.clu, eventTimes, analysis_window, trialGroups);

% figure(5);
% plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);


%% extract LFP
triggers_LFP = uint32(trigger_idx(2:end)*(imec.fsLF/imec.fsAP));
preStim_LFP = -0.2;
postStim_LFP = 3 ;

interval_LFP = [preStim_LFP*imec.fsLF, postStim_LFP*imec.fsLF]; % 1000 samples before plus 999 samples after
snippetSet = imec.readAPSnippetSet(triggers_LFP, interval_LFP);

data = permute(snippetSet.data,[3 2 1]);
% iStart =imec.nGoodChannels-10; iEnd = imec.nGoodChannels;
iStart = 150; 
lfpData_avg = double(data(:,:,iStart));
imagesc((lfpData_avg))



iVal = size(data,3);
% iVal = 77;
for i = 1:iVal
    figure(1);
    imagesc(data(:,:,i))
    title(['Channel ' num2str(i)]);
    caxis([min(min(data(:,:,i)))*min_scale_factor max(max(data(:,:,i)))*max_scale_factor]); colorbar;
pause(1);
end

lfp_kmeans_idx = kmeans(lfpData_avg,4);
first = find(lfp_kmeans_idx==1);
second = find(lfp_kmeans_idx==2);
third = find(lfp_kmeans_idx==3);
fourth = find(lfp_kmeans_idx==4);
[length(first) length(second) length(third) length(fourth)]

plot(lfpData_avg(50,:))

figure(1)
plot(snippetSet.data(5,:,6));
figure(2)
plot(data(6,:,5));

imagesc((lfpData_avg))
max(lfpData_avg(:))
min(lfpData_avg(:))



%%
channel = 129;
scale = [0 1];

data_range = range(data(:,:,channel),'all');
data_image = data(:,:,channel);
imagesc(data_image)
title(['Channel ' num2str(channel)]);
% caxis([min(min(data(:,:,channel)))*min_scale_factor max(max(data(:,:,channel)))*max_scale_factor]); colorbar;
caxis([min(data_image(:))+data_range*scale(1) min(data_image(:))+data_range*scale(2)]); colorbar;

%% Spike amplitude quantification
depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);


%% LFP characterization
myKsDir = bin_file_dir;
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

chanMap = readNPY(fullfile(myKsDir, 'channel_map.npy'));
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

% plot LFP power
dispRange = [0 100]; % Hz
marginalChans = [10:50:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);


%%
% 
% for i = 1:length(trigger_idx_ds)
%     psth(i,:) = [N(trigger_idx_ds(i)-preStimD/bin_size:trigger_idx_ds(i)-1) N(trigger_idx_ds(i):trigger_idx_ds(i)+stimD/bin_size)];
% end
% 
% 
% figure(1)
% bar(time,mean(psth)); ylabel('Spike rate (Hz)'); xlabel('Time relative to stim onset (sec)'); title('All clusters');

%%
% idxWindow = [1:1:30000];
% [data_partial] = imec.readAP_idx(idxWindow);


%%
% myKsDir = 'H:\simon_ephys\Optotest1\_spikeglx_ephysData_probe00_g0';
% nChansInFile = 385;
% syncChanIndex = 385;
% syncDat = extractSyncChannel(myKsDir, nChansInFile, syncChanIndex);
% 
% lfpFs = 2500;  % neuropixels phase3a
% 
% eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);
% 
% sync_b = de2bi(syncDat);
% 

