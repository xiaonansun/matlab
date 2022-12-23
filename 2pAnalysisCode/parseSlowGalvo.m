function [frameStarts, incompleteFrames] = parseSlowGalvo(slowGalvo)

% Identify when frame starts were relative to the analog channels, by
% parsing the slow galvo channel to find oscillatory turnarounds.
% 
% INPUT:
% slowGalvo (vector): voltage trace of the analog MOM channel.

% OUTPUTS:
% frameStarts (vector): index the voltage trace when each frame begins (at the oscillatory turn)
% incompleteFrames (vector): index of the incomplete frames.
% 
% Whenever the recording is stopped, an incomplete frame may or may not be
% recorded. These incomplete frames are given as the second output.
% 
% Goddammit, this is almost the same as samplesPerFrameFromSlowGalvo(). Not
% sure which is better.


%% Parameters

% Minimum change in slow galvo position we'll consider a normal,
% within-frame movement
% minGalvoDiff = 0.15;  % V / ms
% Changed to make this relative, because when we zoom in the range is
% smaller
minGalvoDiff = 0.1 * range(slowGalvo) / (1000 / 30.9);


%% Extract the part of the signal where the voltage is changing
% (clipping off any static portions at the beginning or end)

changes = slowGalvo(1:end-1) ~= slowGalvo(2:end);
firstReal = find(changes, 1);
lastReal = find(changes, 1, 'last');
slowGalvo = slowGalvo(firstReal:lastReal);


%% Find possible frame starts

% Potential frame starts are when either the derivative of the galvo trace
% changes sign, or when there's a too-small change in value (indicating the
% galvo changed direction but was sampled before it passed its previous
% position).
galvoDiff = diff(slowGalvo);
galvoSignChanges = sign(galvoDiff(1:end-1)) ~= sign(galvoDiff(2:end));
% The -1 below is just to make the resulting size of the vector match the
% above, where we have to lose a point
smallGalvoChanges = (abs(galvoDiff(1:end-1)) < minGalvoDiff);
% smallGalvoChanges(end) = 0;
frameStarts = [1 find(galvoSignChanges | smallGalvoChanges)];


%% Remove too-close possible starts

% When our two criteria flagged adjacent samples as possible turnarounds,
% discard the later frame
frameStarts(find(diff(frameStarts) == 1) + 1) = [];


%% Look for incomplete frames due to stopping and restarting

% Find static periods in the slow galvo trace, which indicate that the
% system wasn't actually recording frames
staticPeriods = ~changes(1:end-1) & ~changes(2:end);
staticStarts = find(diff(staticPeriods) == 1) + 1;

% Each time we find a static period, mark the frame before it as incomplete
incompleteFrames = NaN(1, length(staticStarts));
if ~isempty(staticStarts)
  for s = 1:length(staticStarts)
    incompleteFrames(s) = find(sign(frameStarts - staticStarts(s)) < 0, 1, 'last');
  end
end


%% Correct for having clipped front of trace

frameStarts = frameStarts + firstReal - 1;
