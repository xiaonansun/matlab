
% for the source of this script, see: https://github.com/cortex-lab/spikes/blob/master/exampleScript.m
% clear all;
% bin_file_dir = 'J:\simon_ephys\Optotest1\_spikeglx_ephysData_probe00_g1';
% spikeTimes = readNPY([bin_file_dir filesep 'spike_times.npy']);
% spike_clusters = readNPY([bin_file_dir filesep 'spike_clusters.npy']);

%% Loading kilosort data
clear all; close all;

myKsDir = 'J:\simon_ephys\Optotest3\_spikeglx_ephysData_probe00_g0';
myKsDirParts = regexp(myKsDir,filesep,'split');
sp = loadKSdir(myKsDir);

% Computing some useful properties of the spikes and templates
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

%% Drift analysis
descAnalysis = 'Drift Analysis';
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
fDriftAnalysis = figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
fig_configAxis(gca);
strTitle = [myKsDirParts{end} ' ' descAnalysis];
title(strTitle,'Interpreter','none');
saveas(fDriftAnalysis, fullfile(myKsDirParts{1:end-1},[myKsDirParts{end} '_' descAnalysis(~isspace(descAnalysis)) '.pdf']));

%% Quantification of spiking amplitudes
descAnalysis = 'Spiking Amplitudes';
depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
% fSpikingAmplitudes = figure; 
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);
fSpikingAmplitudes = gcf;
fig_configAxis(gca);
strTitle = [myKsDirParts{end} ' ' descAnalysis];
sgtitle(strTitle,'Interpreter','none',...
    'FontSize',12);
saveas(fSpikingAmplitudes, fullfile(myKsDirParts{1:end-1},[myKsDirParts{end} '_' descAnalysis(~isspace(descAnalysis)) '.pdf']));

%% Basic LFP characterization
descAnalysis = 'LFP Characterization';
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

chanMap = readNPY(fullfile(myKsDir, 'channel_map.npy'));
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

dispRange = [0 100]; % Hz
marginalChans = [10:50:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

% fLFP = figure;

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);
fLFP = gcf;
set(fLFP,'Units','inches',...
    'PaperSize',[7 7],...
    'PaperPosition',[0.5 0.5 7 7]);
fig_configAxis(gca);
% strTitle = [myKsDirParts{end} ' ' descAnalysis];
% title(strTitle,'Interpreter','none');
saveas(fLFP, fullfile(myKsDirParts{1:end-1},[myKsDirParts{end} '_' descAnalysis(~isspace(descAnalysis)) '.pdf']));
%% Find opto stim event times

channelMapFile = "C:\Users\Xiaonan Richard Sun\neuropixel-utils\map_files\neuropixPhase3A_kilosortChanMap.mat";

imec = Neuropixel.ImecDataset(myKsDir, 'channelMap', channelMapFile);
sync = imec.readSync();
% mmap = imec.memmapAP_full();
% sync_b = de2bi(sync);

sampling_rate = imec.fsAP;

[trigger_idx,alignedSyncEvents,sync_new] = np_extractSyncEvents(sync, sampling_rate);

eventTimes = trigger_idx/sampling_rate;
alignedSyncEvents = ~alignedSyncEvents;

figure;

% imagesc(alignedSyncEvents);

plot(sum(alignedSyncEvents,1));

alignedSyncEventsBin = np_convertSyncEventsToBinary(alignedSyncEvents);
plot(sum(alignedSyncEventsBin,1));



%% Find trial code categories (four types)
descAnalysis = 'Trial Codes';
endAlignedSyncEvents = align_sync_events_to_end(alignedSyncEvents);

trialType(:,1) = alignedSyncEvents(:,413)==1; %% trials where time point 420 is 1
trialType(:,2) = sum(alignedSyncEvents(:,801:1300),2) ==0; %% trials where there is no signal between 800 and 1300
trialType(:,3) = endAlignedSyncEvents(:,1110) == 1;

trialType(:,4) = ~(trialType(:,1) | trialType(:,2) | trialType(:,3));
% trialType(:,4) = trialTypeRemainder;

% trialType = [trialType1 trialType2 trialType3 trialType4];

fTrialCodes = figure;
imagesc([alignedSyncEvents(trialType(:,1),:);...
    alignedSyncEvents(trialType(:,2),:);...
    alignedSyncEvents(trialType(:,3),:);...
    alignedSyncEvents(trialType(:,4),:)]);
title(descAnalysis)
xlabel('Time sample (30KHz)');
ylabel('Trial')
saveas(fTrialCodes,fullfile(myKsDirParts{1:end-1},[myKsDirParts{end} '_' descAnalysis(~isspace(descAnalysis)) '.pdf']));

%%
% for i = 1:size(alignedSyncEvents,1)
%     Y(i,:) = diff(alignedSyncEvents(i,:));
%     idxChange{i} = find(Y(i,:));
%     idxChangeChange{i} = diff(idxChange{i});
%     idxFirst(i) = find(Y(i,:),1,'first');
%     idxLast(i) = find(Y(i,:),1,'last');
% end


%% Looking at PSTHs aligned to some event

window = [-0.3 4]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
trialGroups = ones(size(eventTimes)); 
trialGroups(trialType(:,2))=2; 
trialGroups(trialType(:,3))=3;
trialGroups(trialType(:,4))=4;
eventTimesByGroup = cell(max(trialGroups),1);
for i = 1:length(eventTimesByGroup)
    eventTimesByGroup{i}=eventTimes(trialType(:,i));
end

% figure;
psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);

%% Divide trials into groups and plot PSTH (i.e. by stimulation parameters)

close all;
clear timeBins depthBins allP normVals meanZScore depthIdx

window = [-1 4]; % look at spike times from 0.3 sec before each event to 1 sec after

depthSurface = 2600; depthTarget = 2200; % target depth in microns
depthDistance = depthSurface - depthTarget;
depthBinSize = 40; % in units of the channel coordinates, in this case µm (default is 80)
ZScoreCutOff = 5;
timeBinSize = 0.01; % seconds
filterDuration = 0.05; % duration of smoothing filter applied
filterLength = filterDuration/timeBinSize;
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'Stimulus onset'; % for figure labeling
eventTimes = eventTimesByGroup; % comment out this line to plot all trials. Otherwise will plot by trial group

if iscell(eventTimes) && numel(eventTimes) > 1
%     meanZScore = nan(length(eventTimes),size(allP,2));
    for i = 1:length(eventTimes)
        [timeBins(i,:), depthBins(i,:), allP{i}, normVals{i}] = psthByDepth(spikeTimes, spikeDepths, ...
            depthBinSize, timeBinSize, eventTimes{i}, window, bslWin);
        depthIdx(i,:) = (depthBins(i,:) >= (depthSurface - depthDistance)) &  (depthBins(i,:) <= depthSurface);
        meanZScore(i,:) = mean(allP{i}(depthIdx(i,:),:));
        stdZScore(i,:) = std(allP{i}(depthIdx(i,:),:));
        meanZScore(i,meanZScore(i,:)>ZScoreCutOff) = nan;
        stdZScore(i,stdZScore(i,:)>ZScoreCutOff) = nan;
        stdZScore(i,:) = stdZScore(i,:)/sqrt(numel((depthSurface-depthDistance):depthBinSize:depthSurface));
        meanZScoreFilt(i,:) = smoothdata(meanZScore(i,:),2,'gaussian',filterLength);
    end
elseif iscell(eventTimes) && numel(eventTimes) == 1
    eventTimes = eventTimes{1};
    [timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
        depthBinSize, timeBinSize, eventTimes, window, bslWin);
else
    [timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
        depthBinSize, timeBinSize, eventTimes, window, bslWin);
end

% Plot mean PSTH
stimDurations = [0.5 1.5];
newXLim = [-0.5 4];
fPSTH = figure(11); hold on;
fPSTH.PaperUnits = 'inches'; fPSTH.PaperSize = [4 3];
fPSTH.Position = [500 500 400 200];
colors = cbrewer('qual','Set1',4);

optoPlotOffset = 1;
bl = boundedline(timeBins(1,:),meanZScoreFilt(1,:),stdZScore(1,:),...
    timeBins(2,:),meanZScoreFilt(2,:),stdZScore(2,:),...
    timeBins(3,:),meanZScoreFilt(3,:),stdZScore(3,:),...
    timeBins(4,:),meanZScoreFilt(4,:),stdZScore(4,:), ...
    'cmap', colors);
for i = 1:length(stimDurations)
    hLine(i) = line([0 stimDurations(i)],[min(cell2mat({bl.YData}))-optoPlotOffset*i min(cell2mat({bl.YData}))-optoPlotOffset*i]);
    hLine(i).Color = 'k';
    hLine(i).LineWidth = 2;
end
for i = 1:length(bl)
    bl(i).LineWidth = 1;
end

strAnimal = 'Optotest';
% myKsDirParts = regexp(myKsDir,filesep,'split');
animal = myKsDirParts{contains(myKsDirParts,strAnimal)}; session = str2num(myKsDirParts{end}(end))+1;
strTitle = {['Opto-inhibition, ' animal ', Session ' num2str(session) ];['1 or 10 mA, Depth: 0 to ' num2str(depthDistance) '\mum']};
xlabel('Time (sec) from stimulus onset');
ylabel('Z-Score')
hTitle = title(strTitle);
hTitle.FontSize = 10;
ax = fig_configAxis(gca);
if exist('newXLim','var'); ax.XLim = newXLim; end

lh = legend(bl);
legnames = {'Cond. 1','Cond. 2', 'Cond. 3', 'Cond. 4'};
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
lh.Location = 'Southeast';
lpos = lh.Position;
lpos(1) = lpos(1) + 0.05;
lh.Position = lpos;

saveas(fPSTH,fullfile(myKsDirParts{1:end-1},[animal myKsDirParts{end} '_PSTH.pdf']))
% print(fPSTH,fullfile(myKsDirParts{1:end-1},[animal myKsDirParts{end} '_PSTH.svg']),'-dsvg')
% print('resize','-');
% hLines = get(gca,'Children');
% legend([hLines(2),hLines(3)],'0.5 sec','1.5 sec',...
%     'Box','off');
%% Plot PSTH at various depth
stimType = {'Cond. 1','Cond. 2','Cond. 3','Cond. 4'};
strAnimal = 'Optotest';
myKsDirParts = regexp(myKsDir,filesep,'split');
animal = myKsDirParts{contains(myKsDirParts,strAnimal)}; session = str2double(myKsDirParts{end}(end))+1;

if ~iscell(allP)
    allPtemp{1} = allP; clear allP
    allP = allPtemp; clear allPtemp;
end

for i = 1:numel(allP)
    %%
    fDepth(i) = figure(5);
    plotPSTHbyDepth(timeBins(i,:), depthBins(i,:), allP{i}, eventName, psthType);
    %         caxis([-2 2]);
    set(fDepth(i),'Units','inches',...
        'PaperSize',[7 7],...
        'PaperPosition',[0.5 0.5 6 6]);
    ax = gca; ax.FontSize = 12;
    ax.YAxis.FontSize=12;
    ax.XAxis.FontSize=12;
    dirParts = regexp(myKsDir,filesep,'split');
    strTitle = {[animal ', Session ' num2str(session) ', ' stimType{i} '. Zero depth is tip of electrode']};
    hTitle(i) = title(strTitle);
    saveas(fDepth(i),fullfile(myKsDirParts{1:end-1},[myKsDirParts{end} '_' num2str(i) '.pdf']));
    print(fDepth(i),fullfile(myKsDirParts{1:end-1},[myKsDirParts{end} '_' num2str(i) '.svg']),'-dsvg');
end

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

%% nested functions

function lineMask = generate_line_mask(se,colLine)
% generates a vertical line spanning 1 or more columns
% colLine = 410:420;
figure;
lineMask = zeros(size(se,1),size(se,2));
lineMask(:,colLine) = ones(size(se,1),length(colLine));
imagesc(se + lineMask);
end

function newSyncEvents = align_sync_events_to_end(se)
%%
flSyncEvents = flip(se,2);
flRealignedSyncEvents = zeros(size(flSyncEvents,1),size(flSyncEvents,2));
% imagesc(flSyncEvents)
[~,iMax] = max(flSyncEvents,[],2);
for i = 1:size(flSyncEvents,1)
flRealignedSyncEvents(i,:)=[flSyncEvents(i,iMax(i):end) flSyncEvents(i,1:iMax(i)-1)];
end
newSyncEvents = flip(flRealignedSyncEvents,2);
% imagesc(flip(flRealignedSyncEvents,2));

end