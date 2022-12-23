function outMat = np_syncReAlignOnTime(inMat, idxTime)
% Realignment of trial code bit signals to a specific point in the trial
% code delivery period
%this is partially complet

%%
idxTime = 30; % comment out this line when executing as afunction

[pos, idx] = max( inMat , [], 2);
pos = logical(pos);
minIdx = min(idx(pos));

% if idxTime 

for i = 1:size(inMat,1)
    if pos(i) == 0
        outMat(i,:) = inMat(i,:);
    else
        if idxTime >= idx(i)
            outMat(i,:) = inMat(i,:);
        else
            outMat(i,:) = [inMat(i,1:idxTime-1) inMat(i,idx(i):end) inMat(i,idxTime:idx(i)-1)];
        end
    end
    
end
imagesc(outMat)
