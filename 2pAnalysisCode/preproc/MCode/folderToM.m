function [MA, MS] = folderToM(folder, stimLen, causal, screenForStability, useTheseModalities, screenForCoinc)
% [MA, MS] = folderToM(folder, stimLen [, causal] [, screenForStability] [, useTheseModalities] [, screenForCoinc])
%
% folder should be the path to a folder containing all_data data structs.
% This function loads each .mat in the folder, converts it to M structs,
% and concatenates over days. The result is two M structs (strength-lumped
% in MA, strength-separated in MS) containing the data over all days from
% the folder.
%
% stimLen determines how long of a stim-locked period to use. Should be
% less than or equal to the shortest imposed wait duration used, or NaN to
% auto-detect for each unit.
%
% causal, if present, determines whether to use a causal filter. Default 0.
%
% screenForStability is default 1.
%
% useTheseModalities, if present and non-empty, limits the results to only
% some modalities (e.g., [-1 0] will yield only auditory and multimodal).
%
% screenForCoinc, if present, determines whether to screen for units with
% too many coincident spikes. 0 means no screening, 1 means to screen all,
% 2 means screen only units with ratings >1. Default 2.


eventFilterSD = 15;


% Optional arguments
if ~exist('causal', 'var')
  causal = 0;
end

if ~exist('screenForStability', 'var')
  screenForStability = 1;
end

if ~exist('useTheseModalities', 'var')
  useTheseModalities = [];
end

if ~exist('screenForCoinc', 'var')
  screenForCoinc = 2;
end

% Find .mat's in the folder
d = dir(folder);

d = d(~[d.isdir]);
d = d(cellfun(@(x) length(x) >= 4 && strcmp(x(end-3:end), '.mat'), {d.name}));


% Loop through files
MA = [];
MS = [];
for i = 1:length(d)
  fprintf('Processing file %02d/%02d, %s   ', i, length(d), d(i).name);
  data = loadOneStruct(fullfile(folder, d(i).name));
  try
    [thisMA, thisMS] = alldataMeans(data, stimLen, NaN, useTheseModalities, eventFilterSD, causal, screenForStability, screenForCoinc);
    
    if isempty(thisMA) || isempty(thisMS)
      fprintf(' No usable units found!\n');
      continue;
    end
    
    % Add filename to M structs
    [thisMA.filename] = deal(d(i).name);
    [thisMS.filename] = deal(d(i).name);
    
    % Concatenate with previous days
    if i == 1
      MA = thisMA;
      MS = thisMS;
    else
      MA = catUnlikeStructs(MA, thisMA);
      MS = catUnlikeStructs(MS, thisMS);
    end
    
    fprintf('\n');
  catch err
    % Print a message if there was a problem loading or processing this file
    fprintf('\nERROR: %s\n', err.message);
    fprintf('**** Dataset not loaded! ****\n');
  end
end


% Print summary info
if length(MA) < 200
  fprintf('%d total units included.\n', length(MA));
else
  % Give yourself a little dopamine
  fprintf('%d total units included. God it feels good to be a gangster.\n', length(MA));
end
fprintf('Done.\n');

