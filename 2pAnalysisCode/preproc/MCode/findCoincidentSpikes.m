function [fracCoincident, expectedFrac, coincident] = findCoincidentSpikes(units, epsilon, totalTimeSpan)
% [fracCoincident, expectedFrac, coincident] = findCoincidentSpikes(units, epsilon, totalTimeSpan)
%
% Compare the fine timing of spike trains form multiple channels to
% determine whether pairs of channels have an excess of nearly coincident
% spikes. Specifically, compare each channel with each other channel, and
% mark a spike that occurs within epsilon ms of a spike on another channel
% as 'coincident'. Also computes the chance level of coincidences,
% approximating each channel as a homogeneous Poisson process (i.e., no
% rate changes, no refractory periods, and continuous instead of discrete
% time).
%
% Inputs:
%   units         -- should be a cell array of length nChannels, with each
%                    cell containing a vector of spike times in ms. Times
%                    should be monotonic, so if using trialized data in
%                    which the clock resets each trial, a greater offset
%                    should be added to each successive trial to prevent
%                    overlapping times.
%   epsilon       -- the jitter to consider spikes 'coincident', in ms.
%                    0.05 is reasonable.
%   totalTimeSpan -- how much time is considered in the spike train given.
%                    E.g., if using trialized data, this should be the sum
%                    of the trial lengths. Should be a scalar.
%
% Outputs:
%   fracCoincident -- vector of length nChannels. Fraction of spikes that
%                     were coincident on this channel.
%   expectedFrac   -- vector of length nChannels. What fraction of spikes
%                     on this channel would occur by chance, if each
%                     channel were a homogeneous Poisson process
%                     rate-matched to the real channel.
%   coincident     -- which spikes coincided with other spikes. Cell array
%                     like units, but each cell contains an array of
%                     logicals of size nUnits x nSpikes, indicating which
%                     unit (if any) each spike coincided with.


%% Handy variables, pre-allocate

nUnits = length(units);

coincident = cell(1, nUnits);
expectedFrac = zeros(1, nUnits);


for u = 1:nUnits
  coincident{u} = false(nUnits, length(units{u}));
end


%% Find coincident spikes

% Loop through units
for u1 = 1:nUnits-1
  
  % Spike train for unit 1. Caching to avoid an extra indexing step.
  s1 = units{u1};
  
  % Loop through all comparison units. Have already compared to units
  % above u1, so only do units below u1.
  for u2 = u1+1:nUnits
    
    % Spike train for unit 2. Caching to avoid an extra indexing step.
    s2 = units{u2};
    
    % Spike indices to compare. Start with first spikes, of course.
    i1 = 1;
    i2 = 1;
    
    % Loop until we run out of spikes on one or the other channel
    while i1 <= length(s1) && i2 <= length(s2)
      if abs(s1(i1) - s2(i2)) < epsilon
        % If spikes are close, mark both as coincident, then move onto the
        % next spike on both channels
        coincident{u1}(u2, i1) = true;
        coincident{u2}(u1, i2) = true;
        
        i1 = i1 + 1;
        i2 = i2 + 1;
      
      elseif s1(i1) < s2(i2)
        % If spikes aren't close, move onto the next spike on the channel
        % that's behind
        i1 = i1 + 1;
        
      else
        % Same, for other channel
        i2 = i2 + 1;
      end
      
    end
    
  end
end



%% Compute fraction coincident

anyCoincident = cell(1, nUnits);
for u = 1:nUnits
  anyCoincident{u} = any(coincident{u});
end
fracCoincident = cellfun(@mean, anyCoincident);


%% Compute chance level of coincidences
% To calculate expectedFrac, use a very quick and dirty algorithm. Assume
% two homogeneous Poisson process (i.e., no rate changes, no refractory
% period). Given the rate of each other unit, compute the probability of
% getting a spike in a given bin that's 2*epsilon wide. Then just multiply
% by the number of spikes (possible windows) in unit 1. This gives the
% expected number of collisions with unit 1. Then, we'll approximate the
% fraction as continuous, to make the last step quick too. Do account for
% different channels having collisions with the same spikes.

% Compute the probability of a coincident spike in a given window,
% 2*epsilon wide, for each channel.
nSpikesPerCh = cellfun(@length, units);
probPerChannel = nSpikesPerCh * 2 * epsilon ./ totalTimeSpan;

% Since coincidences with one channel will fill up a fraction of the
% available spikes, we need to multiply the fractions of slots free, then
% subtract from 1
for u = 1:nUnits
  probPerOtherCh = probPerChannel;
  probPerOtherCh(u) = [];             % all channels except this one
  expectedFrac(u) = 1 - prod(1 - probPerOtherCh);
end



