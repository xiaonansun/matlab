function [alignIdx, trialIdx, frames, baseLength, postLength] = Widefield_getRealignment(fullR, idx, recIdx, trialIdx, recLabels, frames)

%% realign data so baseline is aligned to handle grab and poststim to stimulus
stimOn = sum(fullR(:,ismember(recIdx(~idx),find(ismember(recLabels,{'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'})))),2); %index for stimulus onset in all trials
stimOn = find([0;diff(stimOn)] > 0.5) - 1;

% index for baseline (time before first possible stimulus onset)
baseLength = min(unique(rem(stimOn,frames)))-1;
baseIdx = repmat((0:frames:size(fullR,1)-1)',1,baseLength);
baseIdx = bsxfun(@plus,baseIdx, 1 : baseLength);
baseIdx = baseIdx(:);

% index for post stimulus time
postLength = frames - max(unique(rem(stimOn,frames))); %shortest possible poststim duration
stimIdx = repmat(stimOn, 1, postLength);
stimIdx = bsxfun(@plus,stimIdx, 0 : postLength-1);
stimIdx = stimIdx(:);
alignIdx = sort([baseIdx;stimIdx]);

trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
frames = postLength + baseLength; %new single trial duration in frames
