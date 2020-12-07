function [MA, MS] = identifyGoodUnitsFromM(MA, MS, condsToUse, ratingThresh, minTrCount, SNRThresh)
% [MA, MS] = identifyGoodUnitsFromM(MA, MS [, condsToUse] [, ratingThresh] [, minTrCount] [, SNRThresh])
%
% Keep only "good" units from M structures. These units meet the following
% criteria:
%
% - They have an isolation rating of at least ratingThresh (default 3).
% - They have at least minTrCount trials per condition (default 5).
% - They have an SNR value of at least SNRThresh (default 3.3).
%
% To assess only a subset of conditions for purposes of assessing trial
% counts, specify condsToUse (indexing from 1). For instance, to use only
% multisensory and auditory, which are the second and third conditions, you
% would specify condsToUse as [2 3]. To use all conditions, either do not
% specify condsToUse, or assign it NaN or Inf.

%% Optional arguments

if ~exist('condsToUse', 'var') || isnan(condsToUse) || isinf(condsToUse)
  condsToUse = 1:size(MA(1).cond, 2);
end

if ~exist('ratingThresh', 'var')
  ratingThresh = 3;
end

if ~exist('minTrCount', 'var')
  minTrCount = 5;
end

if ~exist('SNRThresh', 'var')
  SNRThresh = 3.3;
end


%% Do the work

MA_orig = MA;
MS_orig = MS;

nTrials = arrayfun(@(x) min([x.cond(:, condsToUse).nTrials]), MA_orig);

MA = MA_orig([MA_orig.SNR] >= SNRThresh & [MA_orig.rating] >= ratingThresh & nTrials >= minTrCount);
MS = MS_orig([MA_orig.SNR] >= SNRThresh & [MA_orig.rating] >= ratingThresh & nTrials >= minTrCount);
