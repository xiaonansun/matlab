function aligned_sync_events_bin = np_convertSyncEventsToBinary(alignedSyncEvents)

fsPerBit = 6;
aligned_sync_events_bin = zeros(size(alignedSyncEvents,1),round(size(alignedSyncEvents,2)/fsPerBit));
for j = 1:size(alignedSyncEvents,1)
    x = [0 alignedSyncEvents(j,:)];
    idxD = find(diff(x));
    num_of_bins = round(diff(find(diff(x)))/fsPerBit);
    bin_vec = [];
    for i = 1:length(num_of_bins)
        bin_vec = [bin_vec x(idxD(i+1))*ones(1,num_of_bins(i))];
    end
    aligned_sync_events_bin(j,1:length(bin_vec)) = bin_vec;
end

figure;
imagesc(aligned_sync_events_bin);
xlabel('Time sample from first bit (5000 Hz)');
ylabel('Trial #')
title('Aligned trial codes in bit sequence')