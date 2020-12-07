function [smoothedIndex, meanFRs, trialIDs] = instabilityViaMean(alldata, smoothSD, units, maxCut)
% [smoothedIndex, meanFRs, trialIDs] = instabilityViaMean(alldata, smoothSD, units, maxCut)
%
% Compute an index of how unstable the firing rate of each unit is over
% time. Only successful trials are considered. Only the units in 'units'
% are used. 'units' can be given as [] to assess all units.
%
% To compute the index, treat the baseline (-300 to 0 ms from stim onset),
% stimulus (100 to 400 ms after stim onset), and movement epochs (0 to 300
% ms after movement onset) separately. For each, compute the mean for the
% middle portion of the data (from maxCut to 1-maxCut, where maxCut is a
% fraction, e.g., 0.25). The rates for each epoch are then smoothed over
% trials (Gaussian of std dev smoothSD). For each epoch, then, the distance
% of each trial from the 'mid-mean' is taken, and normalized by the
% mid-mean. That is, we get a fractional difference from the mean (of the
% center portion of the data) for each trial for each epoch. The three
% epochs are then averaged for each trial.
%
% smoothedIndex is nUnits x nTrials, and contains the index. We do not
% consider trials that were not successful, did not have neural data, had
% movement durations > 1000 ms, or in which the animal withdrew early.
% These trials get NaN for their index.
%
% meanFRs contains the mean for the three epochs (without subtraction or
% normalization), for considered trials only. This matrix is useful for
% finding trends over the day. The trialId's of considered trials are given
% in trialIDs.


%% Parameters

maxMoveDurFast = 1000;  % ms, for use with no-barrier rats
maxMoveDurSlow = 2000;  % ms, for use with barrier rats


tp.baseStart = -300;
tp.baseEnd = 0;

tp.stimStart = 100;
tp.stimEnd = 400;

tp.moveStart = 0;
tp.moveEnd = 300;



%% Figure out which movement duration cutoff to use

med = nanmedian([alldata.movementDuration]);
if med < 1000
  maxMoveDur = maxMoveDurFast;
else
  maxMoveDur = maxMoveDurSlow;
end



%% Get successful trials

origNTrials = length(alldata);

successes = [alldata.hasNeuralData] == 1 & [alldata.earlyWithdrawal] == 0 & ...
  [alldata.success] == 1 & [alldata.movementDuration] <= maxMoveDur;
alldata = alldata(successes);


%% Determine which units to use

if ~exist('units', 'var') || isempty(units)
  nUnits = length(alldata(1).units);
  units = 1:nUnits;
else
  nUnits = length(units);
end


%% Loop through units, compute the index for each

nTrials = length(alldata);

smoothedIndexGood = NaN(nUnits, nTrials);
meanFRsGood = NaN(nUnits, nTrials);

% Do the real work below!
for u = 1:nUnits
  [smoothedIndexGood(u, :), meanFRsGood(u, :)] = ...
    computeInstabilityIndex(alldata, units(u), tp, smoothSD, maxCut);
end


%% Convert back to original indexing

smoothedIndex = NaN(nUnits, origNTrials);
smoothedIndex(:, successes) = smoothedIndexGood;

meanFRs = NaN(nUnits, origNTrials);
meanFRs(:, successes) = meanFRsGood;

trialIDs = [alldata.trialId];




function [indsSmooth, meanFR] = computeInstabilityIndex(data, unit, tp, smoothSD, maxCut)
% This function does all the real work.

% Find the portion of the data to use for getting the mean
quarter = floor(length(data) * maxCut);
quarters = [quarter length(data)-quarter];

% Handle degenerate case of very short data
if quarter == 0
  indsSmooth = ones(1, length(data));
  return;
end

% Pre-allocate
base = NaN(1, length(data));
stim = NaN(1, length(data));
move = NaN(1, length(data));

% Get counts
for tr = 1:length(data)
  spikes = data(tr).units(unit).spikes;
  
  S = 0;
  M = data(tr).timeInCenter;
  
  base(tr) = sum(spikes > S + tp.baseStart & spikes <= S + tp.baseEnd);
  stim(tr) = sum(spikes > S + tp.stimStart & spikes <= S + tp.stimEnd);
  move(tr) = sum(spikes > M + tp.moveStart & spikes <= M + tp.moveEnd);
end

% To spikes/s, for convenience
base = 1000 * base / (tp.baseEnd - tp.baseStart);
stim = 1000 * stim / (tp.stimEnd - tp.stimStart);
move = 1000 * move / (tp.moveEnd - tp.moveStart);

% Take middle
baseMid = mean(base(quarters(1):quarters(2)));
stimMid = mean(stim(quarters(1):quarters(2)));
moveMid = mean(move(quarters(1):quarters(2)));

% Smooth
baseS = FilterSpikes(smoothSD, base)';
stimS = FilterSpikes(smoothSD, stim)';
moveS = FilterSpikes(smoothSD, move)';

% Compute fractional distance from mean for each epoch
baseIndsS = abs(baseS - baseMid) / (baseMid + 0.01);
stimIndsS = abs(stimS - stimMid) / (stimMid + 0.01);
moveIndsS = abs(moveS - moveMid) / (moveMid + 0.01);

% Compute index
indsSmooth = mean([baseIndsS; stimIndsS; moveIndsS]);
% indsSmooth = FilterSpikes(smoothSD, indsSmooth);


% Output mean rates too
meanFR = mean([base; stim; move]);

% inds = (abs(base - baseLate) + abs(stim - stimLate) + abs(move - moveLate)) ./ ...
%   (baseLate + stimLate + moveLate);


      