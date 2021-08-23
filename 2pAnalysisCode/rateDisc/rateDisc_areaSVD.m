function rateDisc_areaSVD(animal, nrBlocks)

%% some variables
if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
    tPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
    tPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end

if ~exist('nrBlocks', 'var') || isempty(nrBlocks)
    nrBlocks = 49;
end
cPath = [cPath animal filesep];
tPath = [tPath animal filesep 'blockData' filesep];
bOpts.blockDims = 100; %number of dimensions from SVD per block
bOpts.maxLag = 5; %lag for autocorrelation
bOpts.autoConfidence = 0.99; %confidence for autocorrelation test
bOpts.autoThresh = 1.5; %threshold for autocorrelation test
bOpts.snrThresh = 1.6; %threshold for SNR test
bOpts.svdMethod = 'randomized'; %method for blockwise SVD. either 'vanilla' or 'randomized'.
bOpts.nrBlocks = nrBlocks; %keep number of blocks in structure for future reference
bOpts.overlap = 10; %keep number overlap between blocks
bOpts.memLimit = 50; %maximal size per block in gb
bOpts.verbosity = false;
bOpts.trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training

%% get list of recordings and reference image
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
recs = dir([cPath 'SpatialDisc' filesep]);
fprintf('Basepath: %s; Found %d recordings\n', cPath, length(recs));
recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, bOpts.trainingRange); %this sorts recordings by date

fPath = [cPath 'SpatialDisc' filesep recs(1).name filesep]; %Widefield data path
load([fPath filesep 'blueAvg.mat'], 'blueAvg');
load([fPath 'opts2.mat']);
blueAvg = alignAllenTransIm(single(blueAvg),opts.transParams); %align to allen
if opts.frameRate == 15
    fName = 'Vc.mat'; %use original vc
else
    fName = 'rsVc.mat'; %use resampled vc
end

% isolate inner range of allenMask
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask; %mask that is used for all datasets
blueAvg(allenMask) = NaN;
[xRange, yRange] = rateDisc_maskRange(allenMask); % get inner range of allenMask
blueAvg = blueAvg(yRange, xRange);

% get trialcount over all recordings
bhvTrials = cell(1,length(recs));
for iRecs = 1 : length(recs)
    fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
    load([fPath fName], 'bTrials');
    bhvTrials{iRecs} = bTrials;
end

% get index for individual blocks
blockInd = rateDisc_blockInd(blueAvg, nrBlocks, bOpts.overlap);
nrBlocks = length(blockInd); %some blocks might be rejected if outside the mask
save([cPath 'blockInd.mat'],'blockInd');

% estimate size of largest block and determine nr of used frames per recording
nrTrials = size(cat(2,bhvTrials{:}),2); %total nr of trials
nrFrames = (opts.preStim + opts.postStim) * opts.frameRate; %frames per trial
exptSize = size(blockInd{1},1) * 4 * nrTrials * nrFrames / 2^30; %expected size of complete data set in gb.
frameFrac = ceil(exptSize / bOpts.memLimit); %only use fraction of frames to keep memory usage under control
framesPerTrial = floor(nrFrames/frameFrac); %number of frames per trial

% give some feedback
fprintf('First recording: %s; %d recordings - %d trials total \n', fPath, length(recs), nrTrials);
fprintf('MemLimit: %d gb ; Expected size per block: %f gb ; Selected fraction: %d\n', bOpts.memLimit,exptSize,frameFrac);

%% loop through blocks, get data from all recordings and run svd
bU = cell(1, nrBlocks);
if exist([cPath 'blockCnt.txt'],'file')
    delete([cPath 'blockCnt.txt']);
end
for iBlocks = 1 : nrBlocks
    blockData = NaN(size(blockInd{iBlocks}, 1), nrTrials*framesPerTrial, 'single');
    Cnt = 0;
    for iRecs = 1 : length(recs)
        % get data and cut to size
        fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        load([fPath fName], 'Vc', 'U');
        load([fPath 'opts2.mat']);
        U = alignAllenTransIm(single(U),opts.transParams); %align to allen
        U = U(yRange, xRange, :); %crop mask range
        U = reshape(U, [], size(U,3)); %merge pixels
        U = U(blockInd{iBlocks}, :); %isolate block
        
        % subselect frames in each trial
        for iTrials = 1 : size(Vc,3)
            cIdx = find(~isnan(Vc(1,:,iTrials)));
            cIdx = cIdx(randperm(length(cIdx), min([length(cIdx) framesPerTrial])));
            blockData(:, Cnt + (1:length(cIdx))) = U * Vc(:, cIdx, iTrials);
            Cnt = Cnt + length(cIdx);
        end
        clear U Vc
    end
    blockData(:,Cnt+1:end,:) = []; %remove unused entries from blockData
    blockInd{iBlocks}(isnan(mean(blockData,2))) = []; %remove pixels that contain NaNs in any recording
    blockData(isnan(mean(blockData,2)), :) = [];
    bU{iBlocks} = Widefield_compressSVD(blockData, bOpts, false); %compute spatial component
    fID = fopen([cPath 'blockCnt.txt'],'a'); %give some feedback on progress
    fprintf(fID,'Current block: %d/%d; %s \r\n',iBlocks,nrBlocks,datetime); fclose(fID);
end

%% apply U to raw data compute temporal component
for iRecs = 1 : length(recs)   
    % get data and cut to size
    fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
    load([fPath fName], 'Vc', 'U');
    load([fPath 'opts2.mat']);
    U = alignAllenTransIm(single(U),opts.transParams); %align to allen
    U = U(yRange, xRange, :); %crop mask range
    U = reshape(U, [], size(U,3)); %merge pixels
    
    bV = cell(1, length(bU)); %pre-allocate bV
    for iBlocks = 1 : length(bU)
        % compute V for each trial / block
        for iTrials = 1 : size(Vc,3)
            bV{iBlocks}(:, :, iTrials) = bU{iBlocks}' * (U(blockInd{iBlocks}, :) * Vc(:, :, iTrials)); %isolate block
        end
    end
    save([fPath 'bV.mat'], 'bU', 'bV', 'blockInd', 'nrBlocks', '-v7.3'); %save block data
end
if ~exist(tPath, 'dir')
    mkdir(tPath);
end
save([tPath 'bV_' bOpts.trainingRange '.mat'], 'bU', 'bV', 'blockInd', 'nrBlocks', '-v7.3'); %save block data
save([tPath 'svdOpts_' bOpts.trainingRange '.mat'], 'bOpts'); %save option struct

%% do whole frame SVD subsequently
rateDisc_wholeAreaSVD(animal,bOpts.trainingRange)
