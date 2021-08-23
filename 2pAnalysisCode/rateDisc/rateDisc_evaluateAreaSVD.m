function rateDisc_collectEncodingResults(animal,nDims)
%% some variables
% nDims = [100, 500, 1000]; %number of whole-frame components to test
if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
%     cPath = 'Q:\BpodImager\Animals\'; %data path on the server
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
dPath = [cPath animal filesep 'SpatialDisc' filesep]; % path for raw data
trainingRange = 'allAudio';

%% block info
load([bPath 'trialInfo_' trainingRange '.mat'], 'trialCnt', 'recs');
wVfile = matfile([bPath 'wV_' trainingRange '.mat']); %use this for indexing later
lowWVfile = matfile([bPath 'wVlow_' trainingRange '.mat']); %use this for indexing later
load([bPath 'wV_' trainingRange '.mat'],'blockInd');
load([bPath 'wV_' trainingRange '.mat'],'wU');
load([bPath 'bV_' trainingRange '.mat'],'bU');
load([bPath 'mask_' trainingRange '.mat'],'allenMask','xRange','yRange');
mask = allenMask(yRange,xRange);
frameOffset = 50; %skip the first couple of frames to avoid filter artefact
stepSize = 10000;

%get additional dimensions if needed
if size(wU,2) < max(nDims)
    wU = cat(2,wU,lowWVfile.wUlow(:,1:nDims - size(wU,2))); 
end

%% go through recordings and compute variance and prediction error
mse = NaN(length(nDims),length(recs));
for iRecs = 1 : length(recs)
    tic
    fPath = [dPath filesep recs(iRecs).name filesep]; disp(fPath); %session data path
    if exist([fPath 'rsVc.mat'], 'file')
        load([fPath 'rsVc.mat'],'Vc', 'U'); %use downsampled Vc
        opts.frameRate = 15; %adjust framerate
    else
        load([fPath 'Vc.mat'],'Vc', 'U');
    end
    load([fPath 'opts2.mat'],'opts');

    U = alignAllenTransIm(single(U),opts.transParams); %align to allen
    U = U(yRange, xRange, :); %crop to mask range
    
    cIdx = sum(trialCnt(1:iRecs-1))+1 : sum(trialCnt(1:iRecs)); %trials for current recording
    wV = wVfile.wV(:,:,cIdx); %load global dimensions
    if size(wV,1) < max(nDims)
        wV = cat(1,wV,lowWVfile.wVlow(1:nDims - size(wV,1),:,cIdx)); %get additional dimensions
    end

    Vc = Vc(:, 1: size(wV,2),:); %make sure trial lengths are similar
    Vc = reshape(Vc(:,~isnan(Vc(1,:))),size(Vc,1),[]); %exclude NaN frames
    wV = reshape(wV(:,~isnan(wV(1,:))),size(wV,1),[]); %exclude NaN frames

    wV = bsxfun(@minus, wV, mean(wV,2)); %make sure wV is zero-mean
    Vc = bsxfun(@minus, Vc, mean(Vc,2)); %make sure Vc is zero-mean
    
    for iDims = 1 : length(nDims)
        for iFrames = frameOffset : stepSize : size(wV,2)
            if size(wV,2) > iFrames + stepSize
                cData = rateDisc_blockRebuild(mask,blockInd, bU, wV(1:nDims(iDims),iFrames+1 : iFrames + stepSize), wU(:,1:nDims(iDims)));
            else
                cData = rateDisc_blockRebuild(mask,blockInd, bU, wV(1:nDims(iDims),iFrames+1 : end), wU(:,1:nDims(iDims)));
            end
            
            % first run for current recording
            if ndims(U) == 3
                mask = isnan(cData(:,:,1)) | mask;
                U = arrayShrink(U, mask, 'merge');
                varMap = zeros(size(U,1),1);
                errMap = zeros(size(U,1),1);
            end
            cData = arrayShrink(cData, mask, 'merge'); %merge pixels for reconstructed data
            if size(Vc,2) > iFrames + stepSize
                rData = U * Vc(:, iFrames+1 : iFrames + stepSize); %reconstruct 'raw' data
            else
                rData = U * Vc(:, iFrames+1 : end); %reconstruct 'raw' data
            end
            
            varMap = varMap + sum(rData.^2,2); %sum variance in raw data
            errMap = errMap + sum((rData-cData).^2,2); %sum squared prediction error
            clear rData cData
        end
        
        varMap = varMap ./ (size(wV,2) - frameOffset); %variance for each pixel in raw data
        errMap = errMap ./ (size(wV,2) - frameOffset); %squared prediction error in each pixel
        
        cVar = sum(varMap); %total variance
        errVar = sum(errMap); %total squared error
        save([fPath 'blockVar_'  num2str(nDims(iDims)) '.mat'],'cVar','varMap','errVar','errMap');

        mse(iDims,iRecs) = errVar / cVar; %percent error
        fprintf('Used dims: %d; loss: %f percent.\n', nDims(iDims), mse(iDims,iRecs)*100);
    end
    fprintf('Recording %d/%d finished.\n', iRecs, length(recs));
    toc
end

%% save output
cLabel = [];
for iDims = 1 : length(nDims)
    cLabel = [cLabel '_'  num2str(nDims(iDims))];
end
save([bPath 'mse' cLabel '_' trainingRange '.mat'],'mse');