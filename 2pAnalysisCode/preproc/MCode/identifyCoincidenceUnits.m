function [coincUnits, excess] = identifyCoincidenceUnits(alldata, skipUnits, coincThresh, epsilon, printClustering)
% [coincUnits, excess] = identifyCoincidenceUnits(alldata [, skipUnits] [, coincThresh] [, epsilon] [, printClustering])
%
% Figure out which units have too many coincident spikes with other units.
% That is, look for a unit whose spikes occur at almost exactly the same
% time as those of other units. This is done using findCoincidentSpikes()
% (see help for that function for algorithm description).
%
% Inputs:
%   alldata         -- the alldata struct to check
%   skipUnits       -- if not empty, skip these units (should be a vector
%                      of unit numbers). Default: []
%   coincThresh     -- the threshold for coincident spike fraction that
%                      will trigger calling a unit problematic (Default:
%                      0.1)
%   epsilon         -- the spike jitter window in ms (Default: 0.05)
%   printClustering -- whether to print which pairs of units had too many
%                      spikes in common (Default: 0)
%
% Outputs:
%   coincUnits -- indices of units with more coincident spikes than
%                 coincThresh
%   excess     -- what fraction of spikes was in excess of chance, for each
%                 unit (ranges from 0 to 1, or can be a small negative
%                 fraction by chance).


%% Parameters

nTrialsToUse = 100;


%% Optional arguments

if ~exist('skipUnits', 'var')
  skipUnits = [];
end

if ~exist('coincThresh', 'var')
  coincThresh = 0.1;
end

if ~exist('epsilon', 'var')
  epsilon = 0.05;
end

if ~exist('printClustering', 'var')
  printClustering = 0;
end


%% Use only trials with neural data

ad = alldata([alldata.hasNeuralData] == 1);

if isempty(ad)
  error('identifyCoincidenceUnits:noNeuralData', 'No trials with neural data');
end


%% Pick trials to use

% Find center of day and half window length
midTr = floor(length(ad) / 2);
halfTrWindow = floor(nTrialsToUse / 2);
% Find the window around the center of the day
firstTr = midTr - halfTrWindow + 1;
lastTr = midTr + halfTrWindow;
% Ensure the window fits within the data
firstTr = max([1 firstTr]);
lastTr = min([lastTr length(ad)]);

trials = firstTr:lastTr;


%% Convenience variables, pre-allocate

nUnits = length(ad(trials(1)).units);

units = cell(1, nUnits);


%% Compute timing info

% For each trial, want to add an offset so that times are non-overlapping
% (otherwise they restart at 0 for each trial). In addition, we'll throw in
% a 10-second pad between trials, just to be sure.

trialStartTimes = 1000 * arrayfun(@(x) x.parsedEvents.states.state_0(2), ad(trials));
trialOffsets = cumsum([0 diff(trialStartTimes)+10000]);
totalTime = sum([ad(trials).postTrial] - [ad(trials).preTrial]);


%% Assemble 'units' data struct for findCoincidentSpikes()

unitsToCheck = 1:nUnits;
unitsToCheck(skipUnits) = [];

for u = unitsToCheck
  units{u} = arrayfun(@(a, off) a.units(u).spikes' + off, ad(trials), trialOffsets, ...
    'UniformOutput', false);
  units{u} = [units{u}{:}];
end


%% Find coincident spikes

[fracCoinc, expectedFrac, coincident] = findCoincidentSpikes(units, epsilon, totalTime);


%% Identify units to chuck

% Find how many more spikes were coincident than expected by chance.
% Re-scale so that 0 is chance and 1 is all coincident.
excess = (fracCoinc - expectedFrac) ./ (1 - expectedFrac);

coincUnits = find(excess > coincThresh);


%% If requested, print info about which units coincided

if printClustering
  fprintf('\nUnit 1\tUnit 2\tFrac\tTet\tClust\t# conflicts\n');
  for u = coincUnits
    nSpikes = size(coincident{u}, 2);
    nCoincsPerUnit = sum(coincident{u}, 2);
    offenders = find(nCoincsPerUnit' / nSpikes > coincThresh / 2);
    
    if isempty(offenders)
      fprintf('%02d\tNo individual\n', u);
    else
      for off = offenders
        if off == offenders(1)
          fprintf('%02d\t%02d\t%0.2f\t%d\t%d\t%d\n', ...
            u, off, nCoincsPerUnit(off) / nSpikes, ...
            ad(1).units(u).tetrodeNumber, ad(1).units(u).clusterNumber, ...
            length(offenders));
        else
          fprintf('%02d\t%02d\t%0.2f\n', u, off, nCoincsPerUnit(off) / nSpikes);
        end
      end
    end
    
  end
end
