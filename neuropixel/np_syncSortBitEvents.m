function [matA,idxA] = np_syncSortBitEvents(inMat, idxBit)
% Input is one or more cells of vectors or integers referencing the time
% position(s) of the bit

%% 
x = inMat;
% [sel, c] = max( x ~= 0 , [], 2 );
idxA = zeros(size(x,1), length(idxBit));
matA = [];
for i = 1:length(idxBit)
    idxA(:,i) = mean(x(:,idxBit{i}),2)==1;
    idxA = logical(idxA);
    matA = [matA;x(idxA(:,i),:)];
end
