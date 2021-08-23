function frameOut = rateDisc_blockRebuild(mask,blockInd, bU, Vm, wU, blockSelect)
% code to rebuild frames from blockwise or whole frame SVD data.
% mask is used to infer the size of the output frame. blockInd indciates
% the pixel positions in each block, bU is a cell containing spatial
% components for each block, Vm is either a cell with temporal components
% for each block or whole-frame temporal components that can be used to
% compute block-wise scores. In the later case, wU is needed to reconstruct
% the block-wise data first.
% blockSelect is an optional index that can be used to only recounstrct
% selected blocks. In this case only bU{blockSelect} blocks will be
% reconstructed. If blocks were created with spatial overlap, this might
% give a worse reconstruction at the edges of block that are usually
% averaged over.

if ~exist('blockSelect', 'var') || isempty(blockSelect)
    blockSelect = 1 : length(bU);
end
    
[~, cellSize] = cellfun(@size,bU,'UniformOutput',false);
cellSize = cat(2,cellSize{:}); % get number of components in each block

% return whole frame to block components, if needed
if ~iscell(Vm) % Vm is based on whole frame V
    temp = wU * Vm;
    Vm = mat2cell(temp, cellSize, size(temp,2)); %convert to block format
end

% rebuild frame from blocks
dSize = size(Vm{1});
frameOut = NaN(numel(mask), prod(dSize(2:end)),'single');
for iBlocks = blockSelect
    temp = bU{iBlocks} * Vm{iBlocks};
    frameOut(blockInd{iBlocks}, :) = nanmean(cat(3,frameOut(blockInd{iBlocks},:),temp),3);
end
frameOut = reshape(frameOut, size(mask,1), size(mask,2), []);
