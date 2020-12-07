function [MA, MS] = alldataMeans(alldata, stimLen, filterSD, useTheseModalities, eventFilterSD, causal, screenForStability, screenForCoinc)
% [MA, MS] = alldataMeans(alldata [, stimLen] [, filterSD] [, useTheseModalities] [, eventFilterSD] [, causal] [, screenForStability] [, screenForCoinc])
%
% Take one alldata struct and turn it into M structs. MA is the M struct
% with all stimulus strengths averaged, MS has stimulus strengths
% separated.
%
% Trials included meet the following criteria: they are successful, do not
% have an early withdrawal, have neural data, and do not have movement
% times greater than 1000 ms. Each neuron is then scanned for instability
% with markAlldataInstability(). Only stable trials are included if
% screenForStability is 1 (default 1). Units are then screened for whether
% they have an excess of spikes that are coincident with other units. If
% more than 10% of spikes from a unit are coincident with spikes from
% other units, the unit is discarded. If screenForCoinc is 2 instead of 1,
% units with ratings 1 or less are not screened for coincidences. Default
% 2.
%
% Once averaged over trials, the resulting firing rates are double-aligned
% to stimulus and movement onset, then smoothed with a Gaussian with SD
% filterSD (default 50). If causal is 1 (default 0), causal filtering is
% used. If causal is -1, anti-causal filtering is used.
%
% stimLen determines how long of a stimulus-locked epoch will be used.
% Should be less than or equal to the shortest imposed wait duration, or
% NaN to auto-detect for each unit.
%
% This function assumes that stitchData has supplied spike times relative
% to the initial center poke in. The stimAligned field realigns the data
% such that time 0 is the actual stimulus onset.
%
% In addition, if events were synchronized during multi-modal trials, an
% event-triggered average is computed for first events; second, third, and
% fourth events divided by rate; and events 4+ ('late' events). These are
% filtered with an SD of eventFilterSD (default 15).
%
% In addition, an SNR is computed for each unit. This is computed as: (the
% dynamic range) / (the max SEM for any time for any condition plus a weak
% normalizing factor (0.5)).
%
% If useTheseModalities is not empty, it restricts the data to the
% specified modalities. E.g., [-1 1] would use only auditory and visual.


%% Optional arguments

if ~exist('stimLen', 'var')
  timeParams.stimEnd = 500;
else
  timeParams.stimEnd = stimLen;
end

if ~exist('filterSD', 'var') || isnan(filterSD)
  filterSD = 50;
end

if ~exist('useTheseModalities', 'var')
  useTheseModalities = [];
end

if ~exist('eventFilterSD', 'var') || isempty(eventFilterSD)
  eventFilterSD = 15; % ms, NaN to skip event-triggered averaging, generally use 15?
end

if ~exist('causal', 'var')
  causal = 0;
end

if ~exist('screenForStability', 'var')
  screenForStability = 1;
end

if ~exist('screenForCoinc', 'var')
  screenForCoinc = 2;
end


%% Parameters, constants

maxMoveDurFast = 1000;  % ms, for use with no-barrier animals
maxMoveDurSlow = 2000;  % ms, for use with barrier animals

instabilitySmooth = 10;
instabilityThresh = 0.5;
instabilityMaxCut = 0.25;
instabilityMaxTrend = 0.5;

% Codes, must match subfunction
SPIKES_ERROR = 1;
% NOSPIKES = 2;
UNSTABLE = 3;
TOO_MANY_COINC = 4;

% Epoch definitions
timeParams.timeBase = 10;
timeParams.stimStart = -300;
% timeParams.stimEnd = 500;

timeParams.moveStart = -200;
% timeParams.moveEnd = 800;

moveEndFast = 800;
moveEndSlow = 1000;

timeParams.preEvent = -50;
timeParams.postEvent = 200;

% For laser-triggered PETH, if applicable
timeParams.prePulse = -5;
timeParams.postPulse = 20;

removeOptical = 1;



%% Figure out which movement duration cutoff to use

med = nanmedian([alldata.movementDuration]);
if med < 1000
  maxMoveDur = maxMoveDurFast;
  timeParams.moveEnd = moveEndFast;
else
  maxMoveDur = maxMoveDurSlow;
  timeParams.moveEnd = moveEndSlow;
end


%% Get successful trials

successes = [alldata.hasNeuralData] == 1 & [alldata.earlyWithdrawal] == 0 & [alldata.success] == 1;
alldata = alldata(successes);


%% Eliminate trials with overly long movement times

% Add movement duration to data in the process, to make life easy
% for tr = 1:length(alldata)
%   alldata(tr).movementDur = 1000 * (alldata(tr).parsedEvents.states.reward(1) - ...
%     alldata(tr).parsedEvents.states.tone_playing2(2));
% end

alldata = alldata([alldata.movementDuration] <= maxMoveDur);

nUnits = length(alldata(1).units);


%% If present, eliminate optical stimulation trials

if removeOptical && isfield(alldata, 'showOptical')
  alldataOpto = alldata([alldata.showOptical] == 1);
  alldata = alldata([alldata.showOptical] == 0);
end


%% Identify conds

modeByTrial = [alldata.visualOrAuditory];
modes = [-1 0 1];
modeLabels = {'Auditory', 'Multisensory', 'Visual'};

if ~isempty(useTheseModalities)
  goodModes = ismember(modes, useTheseModalities);
  
  modes = modes(goodModes);
  modeLabels = modeLabels(goodModes);
  
  % Trim down alldata to only include the good modes
  alldata = alldata(ismember(modeByTrial, modes));
end



%% Assess stability

if screenForStability
  alldata = markAlldataInstability(alldata, instabilitySmooth, instabilityThresh, instabilityMaxCut, instabilityMaxTrend);
end

%% Identify stimulus strengths and min duration

if any([alldata.nAuditoryEvents] > 0 & [alldata.nVisualEvents] > 0 & ...
    [alldata.nAuditoryEvents] ~= [alldata.nVisualEvents])
  error('Cannot handle cue conflicts');
end

for tr = 1:length(alldata)
  alldata(tr).stimStrength = max(alldata(tr).nAuditoryEvents, alldata(tr).nVisualEvents);
end
strengths = unique([alldata.stimStrength]);
% Below is because Amanda... sometimes?... uses a different field name
if isfield(alldata, 'imposedWaitDuration')
  timeParams.waitDuration = min([alldata.imposedWaitDuration]);
elseif isfield(alldata, 'stimWaitDuration')
  timeParams.waitDuration = min([alldata.stimWaitDuration]);
else
  error('Can''t figure out wait duration field');
end

% Deal with goddamn inconsistency
if timeParams.waitDuration > 0 && timeParams.waitDuration < 2
  timeParams.waitDuration = 1000 * timeParams.waitDuration;
end


%% Metadata, timing

for u = 1:nUnits  
  M(u).tetrode = alldata(1).units(u).tetrodeNumber;
  M(u).cluster = alldata(1).units(u).clusterNumber;
  M(u).rating = alldata(1).units(u).isolation;
  
  %% Separate left/right and modality
  for lr = 1:2
    for mode = 1:length(modes)
      M(u).cond(lr, mode).modality = modes(mode);
      M(u).cond(lr, mode).modeLabel = modeLabels{mode};
      
      if lr == 1
        M(u).cond(lr, mode).leftRight = 'L';
      else
        M(u).cond(lr, mode).leftRight = 'R';
      end
      
      M(u).cond(lr, mode).filterSD = filterSD;

      
      %% Add timing data
      if ~isnan(timeParams.stimEnd)
        M(u).times.stimTimes = timeParams.stimStart : timeParams.timeBase : timeParams.stimEnd;
      else
        M(u).times.stimTimes = timeParams.stimStart : timeParams.timeBase : timeParams.waitDuration;
      end
      M(u).times.moveTimes = timeParams.moveStart : timeParams.timeBase : timeParams.moveEnd;
      M(u).times.eventTimes = timeParams.preEvent : timeParams.postEvent;
      M(u).times.timeBase = timeParams.timeBase;
    end
  end
end


%% Create two clone structures, one for all strengths, one by strength

MA = M;
MS = M;

% Don't waste memory on an extra structure when MA and MS diverge
clear M;


%% Spiking, trial info

badUnits = zeros(1, nUnits);

alldataUntrimmed = alldata;


%% Check for coincident spikes

if screenForCoinc > 1
  skipUnits = find([MA.rating] <= 1 | isnan([MA.rating]));
else
  skipUnits = [];
end

if screenForCoinc > 0
  coincUnits = identifyCoincidenceUnits(alldata, skipUnits);
  % For testing, can uncomment line below and comment line above
  % coincUnits = identifyCoincidenceUnits(alldata, skipUnits, 0.1, 0.05, 1);
  
  if ~isempty(coincUnits)
    badUnits(coincUnits) = TOO_MANY_COINC;
  end
end


%% Loop through units
for u = 1:nUnits
  
  if badUnits(u) ~= 0
    continue;
  end
  
  %% Check for instability, trim trials
  
  if screenForStability
    alldata = getStableTrialsForUnit(alldataUntrimmed, u);

    % Handle unit getting thrown away
    if isempty(alldata)
      % If not, mark bad, skip
      badUnits(u) = UNSTABLE;
      continue;
    end
   
  else
    alldata = alldataUntrimmed;
  end
  
  tp = timeParams;
  tp.stimEnd = MA(u).times.stimTimes(end);
  
  
  
  %% Extract spikes
  
  for lr = 1:2
    for mode = 1:length(modes)
      % OUTDATED: Determine went_left or went_right using left_is_correct
      % and right_is_correct, because later structures lack went_left and
      % went_right.
      if lr == 1
        data = alldata([alldata.visualOrAuditory] == modes(mode) & strcmp({alldata.correctSideName}, 'L'));
      else
        data = alldata([alldata.visualOrAuditory] == modes(mode) & strcmp({alldata.correctSideName}, 'R'));
      end
      
      %% Lump all stimulus strengths
      
      [MA(u).cond(lr, mode).stimAligned, MA(u).cond(lr, mode).moveAligned, ...
        MA(u).cond(lr, mode).stimEvents, badUnit] = ...
        averageData(data, u, tp, filterSD, eventFilterSD, causal);
      
      badUnits(u) = max([badUnits(u) badUnit]);
      
      MA(u).cond(lr, mode).strength = NaN;
      MA(u).cond(lr, mode).RTs = [data.timeInCenter];
      MA(u).cond(lr, mode).moveDurs = [data.movementDuration];
      MA(u).cond(lr, mode).nTrials = length(data);
      
      % Commented out below, because lacking multisensory is not a reason
      % to skip visual.
%       % Check for problem with stimEvents. If found, skip trying in the
%       % future.
%       if isempty(MA(u).cond(lr, mode).stimEvents)
%         eventFilterSD = NaN;
%       end
      
      
      %% Separate by stimulus strength
      
      for s = 1:length(strengths)
        thisData = data([data.stimStrength] == strengths(s));
        
        [MS(u).cond(lr, mode, s).stimAligned, MS(u).cond(lr, mode, s).moveAligned, junk, badUnit] = ...
          averageData(thisData, u, tp, filterSD, eventFilterSD, causal);
        
        badUnits(u) = max([badUnits(u) badUnit]);
        
        MS(u).cond(lr, mode, s).strength = strengths(s);
        MS(u).cond(lr, mode, s).RTs = [thisData.timeInCenter];
        MS(u).cond(lr, mode, s).moveDurs = [thisData.movementDuration];
        MS(u).cond(lr, mode, s).nTrials = length(thisData);
      end
      
    end
  end
end


%% Mark structures as allStrengths or byStrength

[MA.strengthHandling] = deal('allStrengths');
[MS.strengthHandling] = deal('byStrength');


%% Discard bad units

nBad = sum(badUnits > 0);
if nBad == 1
  fprintf('Discarded 1 bad unit:  ');
elseif nBad > 1
  fprintf('Discarded %d bad units: ', nBad);
end

nError = sum(badUnits == SPIKES_ERROR);
nUnstable = sum(badUnits == UNSTABLE);
nCoinc = sum(badUnits == TOO_MANY_COINC);
if nBad > 0
  fprintf('%d unstable, %d coincident, ', nUnstable, nCoinc);
  if nError > 0
    fprintf('%d EXTRACTION ERROR, ', nError);
  end
end
fprintf('%d survive', length(MA) - nBad);


MA(badUnits > 0) = [];
MS(badUnits > 0) = [];

nUnits = length(MA);


%% Add SNR measure

for u = 1:nUnits  
  MA(u).SNR = computeSNR(MA(u).cond);
  MS(u).SNR = computeSNR(MS(u).cond);
end


%% If there were optogenetics trials, compute the laser pulse-triggered PETH

alldata = alldataUntrimmed;

if isfield(alldata, 'showOptical') && ~isempty(alldataOpto)
  units = 1:length(badUnits);
  units(badUnits > 0) = [];
  
  for u = 1:length(units)
    
    % Get stable trials based on alldata
    if screenForStability
      stableTrials = arrayfun(@(ad) ad.units(units(u)).stable, alldata);
      if ~any(stableTrials)
        continue;
      end
      firstId = alldata(find(stableTrials, 1)).trialId;
      lastId = alldata(find(stableTrials, 1, 'last')).trialId;
      goodTrs = ([alldataOpto.trialId] >= firstId & [alldataOpto.trialId] <= lastId);
      if isempty(goodTrs)
        continue;
      end
      
      [PETH, PETHEdges] = computeLaserPETH(alldataOpto(goodTrs), units(u), timeParams);
    else
      [PETH, PETHEdges] = computeLaserPETH(alldataOpto, units(u), timeParams);
    end
      
    MA(u).times.laserPETHEdges = PETHEdges;
    MA(u).laserPETH = PETH;
  end
end


%% Distribute metadata to all strengths for MS

fieldsToCopy = {'modality', 'modeLabel', 'leftRight', 'filterSD'};

for u = 1:nUnits
  for lr = 1:2
    for mode = 1:length(modes)
      nStrengths = size(MS(u).cond, 3);
      if nStrengths > 1
        for strength = 2:nStrengths
          for f = 1:length(fieldsToCopy)
            MS(u).cond(lr, mode, strength).(fieldsToCopy{f}) = MS(u).cond(lr, mode, 1).(fieldsToCopy{f});
          end
        end
      end
    end
  end
end



function [stimAligned, moveAligned, stimEvents, badUnit] = ...
  averageData(data, unit, tp, filterSD, eventFilterSD, causal)

% Codes
SPIKES_ERROR = 1;
% NOSPIKES = 2;


buffer = 100;
nIndivEvents = 4;
lateEvent = 4;
% nConsecutiveEmptyToCauseAbort = 5;


base = tp.timeBase;

bufferPts = buffer / base;
filterPts = filterSD / base;

nTrials = length(data);

badUnit = 0;

% Implement auto-length stim epoch
if isnan(tp.stimEnd)
  tp.stimEnd = tp.waitDuration;
end


%% Trial-average timecourses

% +1 twice because one is for fencepost problem, one is to accommodate
% extra histc bin
stim = zeros(nTrials, (tp.stimEnd - tp.stimStart + 2 * buffer) / base + 1 + 1);
move = zeros(nTrials, (tp.moveEnd - tp.moveStart + 2 * buffer) / base + 1 + 1);

% Find bin boundaries
stimEdges = tp.stimStart - buffer : base : tp.stimEnd + buffer + base;
moveEdges = tp.moveStart - buffer : base : tp.moveEnd + buffer + base;

% Align, bin, filter spikes.
% Also check for many consecutive empty trials. If found, abort -- this
% unit probably wasn't sorted for a while.
for tr = 1:nTrials
  S = data(tr).preStimDelay;
  M = data(tr).timeInCenter;
  
  % It's possible for tetrodes/clusters to have been discarded in the data,
  % but be marked usable in use_these. If we can't access the data, that's
  % what probably happened.
  try
    spikes = data(tr).units(unit).spikes;
  catch
    badUnit = SPIKES_ERROR;
    stimAligned = [];
    moveAligned = [];
    stimEvents = [];
    return;
  end
  
  if ~isempty(spikes)
    sSpikes = ceil(spikes - S);
    stim(tr, :) = histc(sSpikes, stimEdges);
    stim(tr, :) = 1000 / base * FilterSpikes(filterPts, stim(tr, :), causal)';
    
    mSpikes = ceil(spikes - M);
    move(tr, :) = histc(mSpikes, moveEdges);
    move(tr, :) = 1000 / base * FilterSpikes(filterPts, move(tr, :), causal)';
  end
end


% Use non-buffer portion. -1 on tail because of extra histc bin
stimAligned.meanFR = mean(stim(:, bufferPts+1 : end-bufferPts-1), 1);
moveAligned.meanFR = mean(move(:, bufferPts+1 : end-bufferPts-1), 1);

% stdErr. If only one trial, let stdErr = mean
if nTrials > 1
  stimAligned.stdErr = std(stim(:, bufferPts+1 : end-bufferPts-1), 0, 1) / sqrt(nTrials);
  moveAligned.stdErr = std(move(:, bufferPts+1 : end-bufferPts-1), 0, 1) / sqrt(nTrials);
else
  stimAligned.stdErr = stimAligned.meanFR;
  moveAligned.stdErr = moveAligned.meanFR;
end

% Save stim length, enforced wait (use shortest)
stimAligned.shortestStim = min([data.stimDuration]);
stimAligned.shortestImposedWait = tp.waitDuration;


%% Event-average

stimEvents = [];

if ~isnan(eventFilterSD)
  % Use 1 ms timebase
  eventWindow = tp.preEvent:tp.postEvent;
  
  %% Get easy trials with either no auditory or no visual or same aud/vis
  
  usable = arrayfun(@(d) isempty(d.auditoryIsis) || isempty(d.visualIsis) || isequal(d.auditoryIsis, d.visualIsis), data);
  easy = arrayfun(@(d) ~isempty(d.auditoryIsis) && all(d.auditoryIsis(1) == d.auditoryIsis) || ...
    ~isempty(d.visualIsis) && all(d.visualIsis(1) == d.visualIsis), data);
  data = data(usable & easy);
  nTrials = length(data);
  if nTrials == 0
    return;
  end
  
  
  short = zeros(1, nIndivEvents);
  long = zeros(1, nIndivEvents);
    
  firstEvents = NaN(nTrials, length(eventWindow));
  
  rawEventsShortNumbered = NaN(nTrials, nIndivEvents-1, length(eventWindow));
  rawEventsLongNumbered  = NaN(nTrials, nIndivEvents-1, length(eventWindow));
  
  rawEventsShortLate = cell(1, nTrials);
  rawEventsLongLate  = cell(1, nTrials);
  
  
  for tr = 1:nTrials
    
%     %% Check for same auditory and visual stimuli
%     % If not the same, warn and skip
%     if ~isempty(data(tr).auditoryIsis) && ~isempty(data(tr).visualIsis) && ...
%         ~isequal(data(tr).auditoryIsis, data(tr).visualIsis)
%       warning('averageData:UnsynchedStim', ...
%         'Cannot process unsynched stimuli for event-triggered averaging');
%       stimEvents = [];
%       return;
%     end
    
    
    %% Grab the ISI sequence
    if ~isempty(data(tr).auditoryIsis)
      isiSeq = data(tr).auditoryIsis;
    else
      isiSeq = data(tr).visualIsis;
    end
    
    
    %% Find event times
    intLengths = [data(tr).shortInterval, data(tr).longInterval];
    isiLengths = intLengths(isiSeq);
    isiLengths = isiLengths + data(tr).eventDuration;
    % eventTimes is shifted to work with eventWindow
    eventTimes = [0 cumsum(isiLengths)] - tp.preEvent + 1;
    
    
    %% Grab spikes, bin, filter
    % If there's a problem with .spikes, we shouldn't have gotten down here
    % But, spikes can still be empty on a trial or two
%     spikes = ceil(data(tr).units(unit).spikes);
    spikes = ceil(data(tr).units(unit).spikes - data(tr).preStimDelay);
    
    % The below results in 1 ms binned, filtered stim period. The first
    % point in stim is at time tp.preEvent
    dur = eventTimes(end);
    edges = tp.preEvent - buffer : dur + tp.postEvent + buffer;
    if ~isempty(spikes)
      stim = histc(spikes, edges);
      stim = 1000 * FilterSpikes(eventFilterSD, stim, causal)';
      stim = stim(buffer + 1 : end);
      
    else
      stim = zeros(1, length(edges) - 1);
    end
      
      
    %% Numbered events
    
    % First event. Neither short nor long, since no interval precedes it.
    firstEvents(tr, :) = stim(1:length(eventWindow));
    
    % Numbered events 2 through nIndivEvents
    for e = 2:nIndivEvents
      if isiSeq(e-1) == 1
        % If short
        short(e) = short(e) + 1;
        rawEventsShortNumbered(short(e), e-1, :) = stim(eventTimes(e) + eventWindow);
        
      else
        % If long
        long(e) = long(e) + 1;
        rawEventsLongNumbered(long(e), e-1, :) = stim(eventTimes(e) + eventWindow);
      end
    end
    
    
    %% All late events
    
    lateShortEventInds = lateEvent-1 + find(isiSeq(lateEvent-1:end) == 1);
    lateLongEventInds  = lateEvent-1 + find(isiSeq(lateEvent-1:end) == 2);
    
    % Short
    if isempty(lateShortEventInds)
      rawEventsShortLate{tr} = [];
    else
      rawEventsShortLate{tr} = zeros(length(lateShortEventInds), length(eventWindow));
      for e = 1:length(lateShortEventInds)
        rawEventsShortLate{tr}(e, :) = stim(eventTimes(lateShortEventInds(e)) + eventWindow);
      end
    end
    
    % Long
    if isempty(lateLongEventInds)
      rawEventsLongLate{tr} = [];
    else
      rawEventsLongLate{tr} = zeros(length(lateLongEventInds), length(eventWindow));
      for e = 1:length(lateLongEventInds)
        rawEventsLongLate{tr}(e, :) = stim(eventTimes(lateLongEventInds(e)) + eventWindow);
      end
    end
    
  end
  
  %% Average, pack up results
  
  % First events. Pack into both numberedShort and numberedLong
  firstEvent.meanFR = mean(firstEvents);
  firstEvent.stdErr = std(firstEvents) / sqrt(nTrials);
  firstEvent.nEvents = nTrials;
  
  stimEvents.numberedShort(1) = firstEvent;
  stimEvents.numberedLong(1)  = firstEvent;
  
  % Other individually numbered events
  for e = 2:nIndivEvents
    stimEvents.numberedShort(e).meanFR = ...
      squeeze(mean(rawEventsShortNumbered(1:short(e), e-1, :)))';
    stimEvents.numberedShort(e).stdErr = ...
      squeeze(std(rawEventsShortNumbered(1:short(e), e-1, :)) / sqrt(short(e)))';
    stimEvents.numberedShort(e).nEvents = short(e);
    
    stimEvents.numberedLong(e).meanFR = ...
      squeeze(mean(rawEventsLongNumbered(1:long(e), e-1, :)))';
    stimEvents.numberedLong(e).stdErr = ...
      squeeze(std(rawEventsLongNumbered(1:long(e), e-1, :)) / sqrt(long(e)))';
    stimEvents.numberedLong(e).nEvents = long(e);
  end
  
  % Late events
  lateShort = vertcat(rawEventsShortLate{:});
  lateLong  = vertcat(rawEventsLongLate{:});
  
  stimEvents.lateShort.meanFR = mean(lateShort);
  stimEvents.lateShort.stdErr = std(lateShort) / sqrt(size(lateShort, 1));
  stimEvents.lateShort.nEvents = size(lateShort, 1);
  
  stimEvents.lateLong.meanFR = mean(lateLong);
  stimEvents.lateLong.stdErr = std(lateLong) / sqrt(size(lateLong, 1));
  stimEvents.lateLong.nEvents = size(lateLong, 1);
  
  % Other details
  stimEvents.filterSD = eventFilterSD;
end




function SNR = computeSNR(cond)

nConds = numel(cond);

maxMean = 0;
minMean = Inf;
maxErr = 0;

for c = 1:nConds
  maxMean = max([maxMean cond(c).stimAligned.meanFR cond(c).moveAligned.meanFR]);
  minMean = min([minMean cond(c).stimAligned.meanFR cond(c).moveAligned.meanFR]);
  maxErr = max([maxErr cond(c).stimAligned.stdErr cond(c).moveAligned.stdErr]);
end

SNR = (maxMean - minMean) / (maxErr + 0.5);  % Prevent divide-by-0, mildly penalize low FR



function [PETH, binEdges] = computeLaserPETH(ad, unit, timeParams)

binEdges = timeParams.prePulse:timeParams.postPulse;
PETH = zeros(1, length(binEdges));
nPulses = 0;

for tr = 1:length(ad)
  % Use pulse start times. Need to make zero the initial poke, to match up
  % with the spike times, and need to put in ms not seconds
  center = ad(tr).parsedEvents.pokes.C(:, 1);
  center = center(end);
  pulseTimes = 1000 * (ad(tr).parsedEvents.states.optical_pulse(:, 1)' -  center);
  spikeTimes = ad(tr).units(unit).spikes;
  
  for s = spikeTimes'
    % First, get the time of each spike relative to each laser pulse. This
    % tells us what bin it should be in. Let histc do the hard work of
    % figuring whether it falls in the range we care about and which bin it
    % shouyld be in.
    theseCounts = histc(s - pulseTimes, binEdges);
    
    PETH = PETH + theseCounts;
  end
  
  nPulses = nPulses + length(pulseTimes);
end

% Make this an average
PETH = PETH / nPulses;
