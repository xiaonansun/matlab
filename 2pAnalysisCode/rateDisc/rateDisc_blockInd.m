function blockInd = rateDisc_blockInd(refImg, nrBlocks, overlap)
% create index for 'nrBlocks' blocks to fit into 'refImg'. 'overlap'
% determines the nr of pixels of block overlap.

indImg = reshape(1:numel(refImg),size(refImg)); %this is an 'image' with the corresponding indices
blockSize = ceil((size(refImg) + repmat(sqrt(nrBlocks) * overlap, 1, 2))/sqrt(nrBlocks)); %size of each block

Cnt = 0;
colSteps = (0 : blockSize(1) - overlap : size(refImg,1)) + 1; %steps for columns
rowSteps = (0 : blockSize(2) - overlap : size(refImg,2)) + 1; %steps for rows
for iRows = 1 : sqrt(nrBlocks)
    for iCols = 1 : sqrt(nrBlocks)
        % get current block and save index as vector
        colInd = colSteps(iCols) : colSteps(iCols) + blockSize(1) - 1; 
        rowInd = rowSteps(iRows) : rowSteps(iRows) + blockSize(2) - 1;
        
        colInd(colInd > size(refImg,1)) = [];
        rowInd(rowInd > size(refImg,2)) = [];
        
        cBlock = indImg(colInd, rowInd);
        if any(~isnan(refImg(cBlock(:))))  %don't use block if reference contains only NaNs
            Cnt = Cnt + 1;
            blockInd{Cnt} = cBlock(:);
        end
    end
end