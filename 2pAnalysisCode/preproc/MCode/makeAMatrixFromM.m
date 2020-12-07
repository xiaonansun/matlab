function [A, nConds, nStrengths, nTimes, stimLen, normFactors] = makeAMatrixFromM(M, normalize, stimOnly, useTheseConds)
% [A, nConds, nStrengths, nTimes, stimLen, normFactors] = makeAMatrixFromM(M)
% [A, nConds, nStrengths, nTimes, stimLen, normFactors] = makeAMatrixFromM(M, normalize)
% [A, nConds, nStrengths, nTimes, stimLen, normFactors] = makeAMatrixFromM(M, normalize, stimOnly)
% [A, nConds, nStrengths, nTimes, stimLen, normFactors] = makeAMatrixFromM(M, normalize, stimOnly, useTheseConds)
%
% Turn an M struct (produced by alldataMeans or folderToM) into an A
% matrix: nTimepoints*2*nStrengths*nConds x nUnits. The 2 is for left vs.
% right. To recover the tensor, use:
% reshape(A, nTimes, 2, nStrengths, nConds, nUnits)
%
% Inputs:
%
% M              -- the M struct
% normalize      -- whether to normalize columns. 0 is no norm, 1 norms by
%                   column range (across times and conditions), nonzero
%                   values other than exactly 1 norm by range + normalize
%                   ("soft" normalization). Default 1.
% stimOnly       -- whether to use only stimAligned. Default 0.
% useTheseConds  -- whether to use only some conditions. Left and right
%                   will both be used for any used condition. E.g.,
%                   useTheseConds of [1 2] would include auditory and
%                   multimodal for both left and right. This is *not*
%                   indexed based on the modality code; using, e.g., [-1 0]
%                   will produce an error.
%
% Outputs:
%
% A              -- the A matrix
% nConds         -- how many conditions (modalities) there are. if
%                   useTheseConds exists and in non-empty, nConds ==
%                   length(useTheseConds).
% nStrengths     -- how many stimulus strengths there were (left and right
%                   counted separately, with nStrengths the max of the two)
% nTimes         -- number of timepoints
% stimLen        -- number of timepoints in the stimAligned portion per
%                   cond
% normFactors    -- the normalization factors used. NaNs if no
%                   normalization

% Optional params
if ~exist('normalize', 'var')
  normalize = 1;
end

if ~exist('stimOnly', 'var')
  stimOnly = 0;
end

if ~exist('useTheseConds', 'var')
  useTheseConds = [];
end


% Find nConds
if isempty(useTheseConds)
  nConds = size(M(1).cond, 2);
  useTheseConds = 1:nConds;
else
  nConds = length(useTheseConds);
end

% Find nStrengths
nStrengths = size(M(1).cond, 3);

% Find nTimes
stimLens = arrayfun(@(m) length(m.times.stimTimes), M);
stimLen = min(stimLens);
moveLens = arrayfun(@(m) length(m.times.moveTimes), M);
moveLen = min(moveLens);



if stimOnly
  nTimes = stimLen;
  % Pre-allocate
  % *2 is for left vs. right
  A = NaN(nTimes * nConds * nStrengths * 2, length(M));
  
  % Build A
  for u = 1:length(M)
    i = 1;
    for c = 1:nConds
      for s = 1:nStrengths
        for lr = 1:2
          chunk = M(u).cond(lr, useTheseConds(c), s);
          A(i : i + nTimes - 1, u) = chunk.stimAligned.meanFR(1:stimLen)';
          i = i + nTimes;
        end
      end
    end
  end
  
else
  nTimes = stimLen + moveLen;
  % Pre-allocate
  % *2 is for left vs. right
  A = NaN(nTimes * nConds * nStrengths * 2, length(M));
  
  % Build A
  for u = 1:length(M)
    i = 1;
    for c = 1:nConds
      for s = 1:nStrengths
        for lr = 1:2
          chunk = M(u).cond(lr, useTheseConds(c), s);
          A(i : i + nTimes - 1, u) = [chunk.stimAligned.meanFR(1:stimLen) chunk.moveAligned.meanFR]';
          i = i + nTimes;
        end
      end
    end
  end
end


% Norm A
if normalize == 1
  normFactors = range(A);
  normFactors(normFactors == 0) = 1;
  A = bsxfun(@rdivide, A, normFactors);
elseif normalize > 0
  normFactors = range(A) + normalize;
  A = bsxfun(@rdivide, A, normFactors);
else
  normFactors = NaN(1, length(M));
end


