function corrMat = np_sortedCorrelationMap(alignedSyncEvents)

% computes the correlation map by computing the similarity between each
% aligned trial code raw trace
% and then plots the sorted map

corrMat = nan(size(alignedSyncEvents,1));

for i = 1:size(alignedSyncEvents,1)
    for j = 1:size(alignedSyncEvents,1)
    corrMat(j,i) = sum(z(i,:) == z(j,:));
    end
end

figure;
imagesc(sort(sort(corrMat,1),2));