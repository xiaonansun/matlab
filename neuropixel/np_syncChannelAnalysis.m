function np_syncChannelAnalysis(myKsDir)

%%
myKsDir = 'J:\simon_ephys\Optotest4\_spikeglx_ephysData_probe00_g1';
channelMapFile = "C:\Users\Xiaonan Richard Sun\neuropixel-utils\map_files\neuropixPhase3A_kilosortChanMap.mat";
imec = Neuropixel.ImecDataset(myKsDir, 'channelMap', channelMapFile);
sync = imec.readSync();

% finds the index of triggers 
sync_b = de2bi(sync);
% sync_b = int2bit(sync,16); % requires communications toolbox

min_dur = 4; % minimum duration in seconds between each trigger event
interval = imec.fsAP*min_dur;
trigger_idx = 0;

for i=1:length(sync_b)
    if (sync_b(i,1) == 0) && (i-trigger_idx(end) > interval)
        trigger_idx = [trigger_idx i]; % this is the index of the start of individual opto events
    end
end

%% the following script attempts to visualize any clustered formed by
% alignedSyncEvents
mat_sync_events = nan(size(alignedSyncEvents,1)); % first create an empty square matrix
logSE = logical(alignedSyncEvents); % converts alignedSyncEvents to logical
for i = 1:size(logSE,1)
    for j = i:size(logSE,1)
    mat_sync_events(j,i) = sum(logSE(j,:) & logSE(i,:)); % computes a matrix with the sum of the products of every combination of vectors
    end
end
% plot(sort(mat_sync_events(:,1)),'.k')
% Unfortunately no clusters were resolved with this analysis

imec.inspectAP_timeWindow([0 240],'channels',385)

diff(unique(sync)')

%% load synchronization data

syncChanIndex = 385;
syncDat = extractSyncChannel(myKsDir, nChansInFile, syncChanIndex);

eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);

% - eventTimes{1} contains the sync events from digital channel 1, as three cells: 
% - eventTimes{1}{1} is the times of all events
% - eventTimes{1}{2} is the times the digital bit went from off to on
% - eventTimes{1}{2} is the times the digital bit went from on to off

% To make a timebase conversion, e.g. between two probes:
% [~,b] = makeCorrection(syncTimesProbe1, syncTimesProbe2, false);

% and to apply it:
% correctedSpikeTimes = applyCorrection(spikeTimesProbe2, b);

%%
% Additional analysis of the sync channel
% The first part of the code attempts multiple permutations of shifting the
% bit-wise codes to decode characters.
% Context: teensy sent bits containing character information to channel
% 385 (the sync channel) of the neuropix phase 3A probe. This is an attempt
% to decode the individual characters. So far, this effort has been
% unsuccessful.

%%
bin_shift = 20;
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