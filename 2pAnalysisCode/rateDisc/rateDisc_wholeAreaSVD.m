function rateDisc_wholeAreaSVD(animal, trainingRange)

%% some variables
if ispc
    cPath = 'U:\smusall\BpodImager\Animals\';
    tPath = 'U:\smusall\BpodImager\Animals\';
%     tPath = 'X:\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
    tPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end

cPath = [cPath animal filesep];
tPath = [tPath animal filesep 'blockData' filesep];
bOpts.blockDims = 1000; %number of dimensions from whole brain SVD
bOpts.maxLag = 5; %lag for autocorrelation
bOpts.autoConfidence = 0.99; %confidence for autocorrelation test
bOpts.autoThresh = 1.5; %threshold for autocorrelation test
bOpts.snrThresh = 1.6; %threshold for SNR test
bOpts.svdMethod = 'randomized'; %method for blockwise SVD. either 'vanilla' or 'randomized'.
bOpts.overlap = 10; %keep number overlap between blocks
bOpts.memLimit = 60; %maximal size per block in gb
bOpts.verbosity = false;
bOpts.trialDur = 90; %expected number of frames per trial
bOpts.trainingRange = trainingRange; %use this to run analysis only in a certain range of training
bOpts.nDims = 500; %number of dimensions in new V/U

%% get list of recordings and reference image
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
recs = dir([cPath 'SpatialDisc' filesep]);
recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, bOpts.trainingRange); %this sorts recordings by date

% get trialcount over all recordings
bhvTrials = cell(1,length(recs));
for iRecs = 1 : length(recs)
    fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
    try
        load([fPath 'opts2'], 'opts'); %check sampling rate
        if opts.frameRate == 15
            load([fPath 'Vc'], 'bTrials'); %check for Vc, dont use recording otherwise
        else
            load([fPath 'rsVc'], 'bTrials'); %check for Vc, dont use recording otherwise
        end
        bhvTrials{iRecs} = bTrials;
        trialCnt(iRecs) = length(bTrials);
        useIdx(iRecs) = true;
    catch
        useIdx(iRecs) = false;
    end
end
bhvTrials = bhvTrials(useIdx);
trialCnt = trialCnt(useIdx);
recs = recs(useIdx);
save([tPath 'trialInfo_' bOpts.trainingRange '.mat'], 'bhvTrials', 'trialCnt','recs'); %save trial info
nrTrials = sum(trialCnt); %total nr of trials

load([tPath 'bV_' bOpts.trainingRange '.mat'], 'nrBlocks'); %get number of blocks
bOpts.nrBlocks = nrBlocks; %keep number of blocks in structure for future reference

%% loop through blocks, get data from all recordings and run svd
allV = cell(1, nrBlocks);
dimCnt = NaN(1,nrBlocks);
for iBlocks = 1 : nrBlocks
    Cnt = 0; 
    for iRecs = 1 : length(recs)
        fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        cFile = matfile([fPath 'bV.mat']); %blockV file
        temp = cFile.bV(1,iBlocks); %get current block from mat file
        
        if iRecs == 1
            allV{1,iBlocks} = NaN(size(temp{1},1),bOpts.trialDur,nrTrials, 'single'); %pre-allocate array for current block
        end
        allV{1,iBlocks}(:,:,Cnt + (1:size(temp{1},3))) = temp{1}(:,1:bOpts.trialDur,:); %collect in cell array
        Cnt = Cnt + size(temp{1},3);
    end
    dimCnt(iBlocks) = size(allV{1,iBlocks},1);
    disp([num2str(iBlocks) '-' num2str(iRecs)]);
end

%% apply U to raw data compute temporal component
allV = cat(1,allV{:}); %combine all blocks
allV = reshape(allV,size(allV,1),[]); %merge all frames
useIdx = ~isnan(allV(1,:)); %find NaN frames
allV = allV(:, useIdx); %reject from data

%run svd over all blocks and separate signal from noise dimensions
[wU, wV1, wUlow, wVlow1] = Widefield_compressSVD(allV, bOpts, true); clear allV

% pre-allocate temporal dims that includes NaN frames
wV = NaN(size(wV1,1), length(useIdx), 'single');
wV(:, useIdx) = wV1; clear wV1
wV = reshape(wV, size(wV,1), [], nrTrials);

wVlow = NaN(size(wVlow1,1), length(useIdx), 'single');
wVlow(:, useIdx) = wVlow1; clear wVlow1
wVlow = reshape(wVlow, size(wVlow,1), [], nrTrials);

% save down global dimensions
if ~exist(tPath, 'dir')
    mkdir(tPath);
end
blockInd = cFile.blockInd;
save([tPath 'wV_' bOpts.trainingRange '.mat'], 'wU', 'wV', 'dimCnt', 'blockInd', 'nrBlocks', '-v7.3'); %save block data
save([tPath 'wVlow_' bOpts.trainingRange '.mat'], 'wUlow', 'wVlow', 'dimCnt', 'blockInd', 'nrBlocks', '-v7.3'); %save block data

% also save some mask info
% isolate inner range of allenMask
load('allenDorsalMapSM.mat', 'dorsalMaps')
allenMask = dorsalMaps.allenMask; %mask that is used for all datasets
[xRange, yRange] = rateDisc_maskRange(allenMask); % get inner range of allenMask
bOpts.nDims = min([bOpts.nDims,size(wV,1)]); %make sure there are enough dims for later analysis
save([tPath 'globalOpts_' bOpts.trainingRange '.mat'], 'bOpts', 'recs', 'bhvTrials'); %save option struct

%% make new U
load([tPath 'bV_' bOpts.trainingRange '.mat'],'bU', 'blockInd');
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
newU = blockU * wU(:,1:bOpts.nDims); %create newU that combines block-wise and whole page components

blockU = arrayShrink(blockU,mask,'split');
newU = arrayShrink(newU,mask,'split');

save([tPath 'blockU_' bOpts.trainingRange '.mat'], 'blockU', 'newU', '-v7.3');
save([tPath 'mask_' bOpts.trainingRange '.mat'],'allenMask','mask','xRange','yRange'); %add mask across sessions here

%% do qr for locaNMF code
wV = wV(1:bOpts.nDims,:,:);
wV = reshape(wV, bOpts.nDims, []); %make sure wV is in 2D
nanIdx = isnan(wV(1,:)); %keep index for NaN frames
wV(:, nanIdx) = []; %remove NaNs
[q, r] = qr(wV', 0); %get qr results
save([tPath 'wQR_' bOpts.trainingRange '.mat'], 'q', 'r', 'nanIdx', '-v7.3');
