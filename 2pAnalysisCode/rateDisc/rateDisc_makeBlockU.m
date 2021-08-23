function rateDisc_makeBlockU(animal, ndims)
% code to combine tiled 'bU' back into a single U matrix that is composed
% of individual blocks. This can be given to the locaNMF code.

%% some variables
if ~exist('ndims', 'var') || isempty(ndims)
    ndims = 500;
end
if ispc
%     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for global dimensions

load([bPath 'wV_allAudio.mat'],'wU');
load([bPath 'bV_allAudio.mat'],'bU', 'blockInd');
load([bPath 'mask.mat'],'allenMask','xRange','yRange');
mask = allenMask(yRange,xRange);

[~, cellSize] = cellfun(@size,bU,'UniformOutput',false);
cellSize = cat(2,cellSize{:}); % get number of components in each block

% rebuild block-wise U from individual blocks
blockU = zeros(numel(mask), sum(cellSize),'single');
edgeNorm = zeros(numel(mask),1,'single');
Cnt = 0;
for iBlocks = 1 : length(bU)
    cIdx = Cnt + (1 : size(bU{iBlocks},2));
    blockU(blockInd{iBlocks}, cIdx) = blockU(blockInd{iBlocks}, cIdx) + bU{iBlocks};
    edgeNorm(blockInd{iBlocks}) = edgeNorm(blockInd{iBlocks}) + 1;
    Cnt = Cnt + size(bU{iBlocks},2);
end
edgeNorm(edgeNorm == 0) = 1; %remove zeros to avoid NaNs in blockU

%normalize blockU by dividing pixels where blocks overlap
blockU = bsxfun(@rdivide, blockU, edgeNorm);
blockU = reshape(blockU, size(mask,1), size(mask,2), []);
mask = mask | sum(blockU,3) == 0; %don't keep zero pixels (these are pixels that are inconsistent across sessions);
blockU = arrayShrink(blockU,mask,'merge');
newU = blockU * wU(:,1:ndims); %create newU that combines block-wise and whole page components

blockU = arrayShrink(blockU,mask,'split');
newU = arrayShrink(newU,mask,'split');

save([bPath 'blockU.mat'], 'blockU', 'newU', '-v7.3');
save([bPath 'mask.mat'],'allenMask','mask','xRange','yRange'); %add mask across sessions here


