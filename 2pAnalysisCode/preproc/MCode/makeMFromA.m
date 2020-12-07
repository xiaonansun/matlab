function newM = makeMFromA(A, M, nConds, nStrengths, nTimes, stimLen, stdErrA)
% newM = makeMFromA(A, M, nConds, nStrengths, nTimes, stimLen)
% newM = makeMFromA(A, M, nConds, nStrengths, nTimes, stimLen, stdErrA)
%
% Make a (new) M structure from an A matrix. M should be an example M
% struct to pull metadata from. If M has the same number of units as A,
% makeMFromA will preserve the metadata for each unit. Otherwise, all
% metadata will come from M(1).
%
% If supplied and the same size as A, stdErrA will be used to fill in
% .stdErr fields in newM. Otherwise, .stdErr will be [].


% Optional argument
if ~exist('stdErrA', 'var') || ~isequal(size(A), size(stdErrA))
  stdErrA = [];
end



nUnits = size(A, 2);

% Prepare newM
if length(M) == nUnits
  % Can preserve unit-wise metadata
  newM = M;
else
  % Use M(1) for all metadata
  newM = repmat(M(1), 1, nUnits);
end


for u = 1:nUnits
  % Reshape A to: nTimes x 2 x S x C
  meanResp = reshape(A(:, u), nTimes, 2, nStrengths, nConds);
  
  
  % Fill in .meanFR fields
  for c = 1:nConds
    for s = 1:nStrengths
      for lr = 1:2
        newM(u).cond(lr, c, s).stimAligned.meanFR = squeeze(meanResp(1:stimLen, lr, s, c))';
        newM(u).cond(lr, c, s).moveAligned.meanFR = squeeze(meanResp(stimLen+1:nTimes, lr, s, c))';
      end
    end
  end
  
  
  % Fill in .stdErr fields
  if ~isempty(stdErrA)
    % Use data from stdErrA
    % Reshape stdErrA to: nTimes x 2 x S x C
    stdErr = reshape(stdErrA(:, u), nTimes, 2, nStrengths, nConds);
    for c = 1:nConds
      for s = 1:nStrengths
        for lr = 1:2
          newM(u).cond(lr, c, s).stimAligned.stdErr = squeeze(stdErr(1:stimLen, lr, s, c))';
          newM(u).cond(lr, c, s).moveAligned.stdErr = squeeze(stdErr(stimLen+1:nTimes, lr, s, c))';
        end
      end
    end
    
  else
    % Use []
    for c = 1:nConds
      for s = 1:nStrengths
        for lr = 1:2
          newM(u).cond(lr, c, s).stimAligned.stdErr = [];
          newM(u).cond(lr, c, s).moveAligned.stdErr = [];
        end
      end
    end
  end
  
end

