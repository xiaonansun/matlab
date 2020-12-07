% clear all;
% bin_file_dir = 'H:\simon_ephys\Optotest1\_spikeglx_ephysData_probe00_g1';
% spikeTimes = readNPY([bin_file_dir filesep 'spike_times.npy']);
% spike_clusters = readNPY([bin_file_dir filesep 'spike_clusters.npy']);

%% Loading kilosort data
clear all; close all;

myKsDir = 'H:\simon_ephys\Optotest4\_spikeglx_ephysData_probe00_g1';
sp = loadKSdir(myKsDir)

spikeTimes = sp.st;

%% Analyzing drift
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);


%% Quantification of spiking amplitudes
depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

figure;
[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);

%% Basic LFP characterization

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

figure;
plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);


%% Computing some useful properties of the spikes and templates
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);


%% Find opto stim event times

channelMapFile = "C:\Users\Xiaonan Richard Sun\neuropixel-utils\map_files\neuropixPhase3A_kilosortChanMap.mat";

imec = Neuropixel.ImecDataset(myKsDir, 'channelMap', channelMapFile);
sync = imec.readSync();
mmap = imec.memmapAP_full();
sync_b = de2bi(sync);

min_dur = 4; % minimum duration in seconds between each trigger event
interval = imec.fsAP*min_dur;
trigger_idx = 0;

for i=1:length(sync_b)
    if (sync_b(i,1) == 0) && (i-trigger_idx(end) > interval)
        trigger_idx = [trigger_idx i];
    end
end

eventIdx = trigger_idx(2:end);
eventTimes = trigger_idx(2:end)/sp.sample_rate;

aligned_sync_window = [0 0.05]; % window for aligning sync events into a matrix
alignedSyncEvents = zeros(length(eventIdx),abs(aligned_sync_window(1)*sp.sample_rate-1)+aligned_sync_window(2)*sp.sample_rate);

for i = 1:length(eventIdx)
   alignedSyncEvents(i,:) = 1-[sync_b(eventIdx(i)+aligned_sync_window(1)*sp.sample_rate:eventIdx(i)-1,1)' sync_b(eventIdx(i):eventIdx(i)+aligned_sync_window(2)*sp.sample_rate,1)'];
end

plot(aligned_sync_window(1)*sp.sample_rate:1:aligned_sync_window(2)*sp.sample_rate, alignedSyncEvents(1,:))
imagesc(alignedSyncEvents);

%%
bin_shift = 0;
sync_b1 = circshift(alignedSyncEvents(1,:),bin_shift);
minibin = 25;
subIdx = [1 7 13 19];
inputBin = [];
for i = 1:floor(length(sync_b1)/minibin)
    temp = sync_b1(1+(i-1)*minibin:i*minibin);
    inputBin = [inputBin temp(subIdx)];
end

shift = 0;
% iB = circshift(flip(inputBin(1:floor(length(inputBin)/8)*8)),shift);

iB = circshift(inputBin(1:floor(length(inputBin)/8)*8),shift);
iBrshp = reshape(iB',[8,length(iB)/8]);
C = bi2de(iBrshp');
char(C)'
%%
for i = 1:size(alignedSyncEvents,1)
    Y(i,:) = diff(alignedSyncEvents(i,:));
    idxChange{i} = find(Y(i,:));
    idxChangeChange{i} = diff(idxChange{i});
    idxFirst(i) = find(Y(i,:),1,'first');
    idxLast(i) = find(Y(i,:),1,'last');
end


% [IDX, C] = kmeans(alignedSyncEvents, 4);
% [length(find(IDX ==1)) length(find(IDX ==2)) length(find(IDX ==3)) length(find(IDX ==4))]
%% Looking at PSTHs aligned to some event

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
trialGroups = ones(size(eventTimes)); 

figure;
psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);

%% plot by depth
close all;

depthBinSize = 80; % in units of the channel coordinates, in this case µm
timeBinSize = 0.01; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, eventTimes, window, bslWin);

hDepth = figure(5);
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);
set(hDepth,'Units','inches',...
    'PaperPosition',[0.5 0.5 3 3]);
ax = gca; ax.FontSize = 8; 
ax.YAxis.FontSize=4;
ax.XAxis.FontSize=4;
dirParts = regexp(myKsDir,'\','split');
print([fullfile(dirParts{1:end-1}) filesep dirParts{end-1} dirParts{end} '.pdf'],'-dpdf','-bestfit');

% saveas(hDepth,[fullfile(dirParts{1:end-1}) filesep dirParts{end-1} dirParts{end} '.pdf']);

%% Loading raw waveforms (or spike-triggered LFP, e.g.)

gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==155);

wf = getWaveForms(gwfparams);

figure; 
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;