function alldata = markAlldataInstability(alldata, smoothSD, thresh, maxCut, maxTrendSize)
% alldata = markAlldataInstability(alldata)
% alldata = markAlldataInstability(alldata, smoothSD)
% alldata = markAlldataInstability(alldata, smoothSD, thresh)
% alldata = markAlldataInstability(alldata, smoothSD, thresh, maxCut)
% alldata = markAlldataInstability(alldata, smoothSD, thresh, maxCut, maxTrendSize)
%
% Assesses stability of firing rate over trials, for each unit in alldata.
% Uses instabilityViaMean() to produce an instability index (II) for each
% trial for each unit -- see 'help instabilityViaMean' for description of
% that algorithm.
%
% Once the index is computed, for each unit, checks if there are regions of
% instability near the beginning or end of the data. If there are at least
% 2*smoothSD points where the II exceeds thresh, and the region lies
% entirely within maxCut of the beginning or end of the data (e.g., 0.25),
% marks the points from the nearer end of the data through the unstable
% region as unstable. If the unstable region is in the middle of the data
% (i.e., more than maxCut from either end), marks all trials as unstable
% for this unit.
%
% If end-clipping is successful or unneeded, checks whether there is an
% overall trend in the data. This is done by regressing a line to the data
% (only in the preserved trials), using the meanFR output value from
% instabilityViaMean(). The values at the ends of the line are computed,
% their ratio is taken (larger / smaller), and 1 is subtracted. Thus, if
% there is no trend in the data, the 'trend size' is 0, while if the firing
% rate rises (or falls) by a factor of 2 (100%), the trend size is 1. If
% the trend size is greater than maxTrendSize, the unit is marked unstabe
% for all trials.
%
% Optional arguments:
%
% smoothSD     -- Passed to instabilityViaMean(), and used for finding
%                 long-enough zones of instability. Default 10
% thresh       -- how high the II must be before we consider it potentially
%                 unstable. Default 0.5
% maxCut       -- The maximum fraction of the data to clip, from either
%                 end. Default 0.25
% maxTrendSize -- The maximum trend size before the whole unit is marked
%                 unstable. Default 0.4
%
%
% Output alldata struct modifications:
%
% The .units field contains two new subfields:
%  .stable is 1 or 0, and indicates whether that unit was stable for this
%      trial.
%  .stabilityVals is a struct, and contains the II (different for each
%      trial for each unit) and trend size (same for all trials for a given
%      unit).
%
% In addition, the parameters used for assessing stability are in a new
% subfield at the top level, called .stabilityParams.

%% Optional arguments

if ~exist('smoothSD', 'var')
  smoothSD = 10;
end

if ~exist('thresh', 'var')
  thresh = 0.5;
end

if ~exist('maxCut', 'var')
  maxCut = 0.25;
end

if ~exist('maxTrendSize', 'var')
  maxTrendSize = 0.5;
end



%% Parameters

% How many points must exceed thresh before we try to clip.
% Value multiplied by smoothSD!
nInARow = 2;


%% Error check

if maxCut >= 0.5
  error('maxCut must be < 0.5');
end


%% Convenience variables

firstGoodTr = find([alldata.hasNeuralData], 1);
nUnits = length(alldata(firstGoodTr).units);


%% Compute assessment for all units

[smoothedIndex, meanFRs, trialIDs] = instabilityViaMean(alldata, smoothSD, [], maxCut);

alldataInds = find(ismember([alldata.trialId], trialIDs));


%% Find violators

nonNans = find(~isnan(smoothedIndex(1, :)));
overThresh = (smoothedIndex(:, nonNans) >= thresh);


%% Prepare to quickly find points to clip

% Boxcar filter each unit
cs = cumsum(overThresh, 2);
meanOverThresh = cs(:, nInARow*smoothSD+1:end) - cs(:, 1:end-nInARow*smoothSD);

% Find indices of maxCut and end-maxCut. Call them 'quarters'
half = floor(size(overThresh, 2) / 2);
quarter = floor(size(overThresh, 2) * maxCut);
quarters = [quarter length(overThresh)-quarter];

meanFRNoNans = meanFRs(:, nonNans);


for u = 1:nUnits
  
  %% Find ends to clip
  
  % Efficiently find when there are nInARow*smoothSD violators in a row
  
  % First half: find boxcar totally full of 1's
  % Note that names are "reversed", because we're looking for problems
  lastProblemNonNans = find(meanOverThresh(u, 1:half) == nInARow*smoothSD, 1, 'last') + nInARow*smoothSD;
  
  % Second half: find boxcar totally full of 1's
  firstProblemNonNans = find(meanOverThresh(u, half+1:end) == nInARow*smoothSD, 1) + half + 1;
  
  
  %% Handle instability, integrate into alldata
  
  % Prep stability marks, instability index, and trendSize
  for tr = 1:length(alldata)
    alldata(tr).units(u).stable = false;
    alldata(tr).units(u).stabilityVals.instabilityIndex = NaN;
    alldata(tr).units(u).stabilityVals.trendSize = NaN;
  end
  
  for tr = alldataInds
    alldata(tr).units(u).stabilityVals.instabilityIndex = smoothedIndex(u, tr);
  end
  
  % Check whether the unit can be salvaged. If so, mark ok trials. If not,
  % leave everything as-is, with no ok trials marked.
  if (isempty(lastProblemNonNans)  || lastProblemNonNans <= quarters(1)) && ...
     (isempty(firstProblemNonNans) || firstProblemNonNans > quarters(2))
    
    % Check whether to clip ends
    if isempty(lastProblemNonNans)
      lastProblemNonNans = 0;
      firstInd = 1;
    else
      firstInd = alldataInds(lastProblemNonNans+1);
    end
    if isempty(firstProblemNonNans)
      firstProblemNonNans = size(meanFRNoNans, 2) + 1;
      lastInd = length(alldata);
    else
      lastInd = alldataInds(firstProblemNonNans-1);
    end
    
    % Check for linear trends in the mean
    meanFRThisUnit = meanFRNoNans(u, lastProblemNonNans+1 : firstProblemNonNans-1);
    p = polyfit(1:length(meanFRThisUnit), meanFRThisUnit, 1);
    endVals = polyval(p, [1 length(meanFRThisUnit)]);
    trendSize = endVals(1) / endVals(2);
    if abs(trendSize) < 1
      trendSize = 1 / trendSize;
    end
    trendSize = trendSize - 1;
    
    % Save trendSize
    for tr = 1:length(alldata)
      alldata(tr).units(u).stabilityVals.trendSize = trendSize;
    end
    
    % Check that trend is small enough, mark ok trials as ok
    if trendSize < maxTrendSize
      for tr = firstInd : lastInd
        alldata(tr).units(u).stable = true;
      end
    end
    
  end
end


% Save parameters for deciding instability
for tr = 1:length(alldata)
  alldata(tr).stabilityParams.instabilityThresh = thresh;
  alldata(tr).stabilityParams.instabilitySmoothingSD = smoothSD;
  alldata(tr).stabilityParams.trendSizeThresh = maxTrendSize;
end

