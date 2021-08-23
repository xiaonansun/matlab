function rateDisc_videoRebuild(cPath,fileExt)
% code to rebuild behavioral video data based on beta coefficients from
% lienar model. The reconstructed video data should show which spatial
% patterns in the behavioral video are used to inform the model about
% changes in widefield data.

if ~exist('fileExt','var')
    fileExt = '';
end
    
if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

undoQR = true;  % default is orthogonalized video. Undo that QR step.
if strcmpi(fileExt,'org') || strcmpi(fileExt,'vidOnly') %cases were unQRed video was used
    undoQR = false;
end

%% load data
load([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']); %load combined SVD
load([cPath fileExt 'regData.mat'],'regIdx', 'regLabels', 'rejIdx','redQRR'); %load R from QR decomposition
load([cPath fileExt 'dimBeta.mat'],'dimBeta'); %load model betas for orthogonalized design matrix

load([cPath 'BehaviorVideo' filesep 'segInd1.mat'],'ind'); %load segment index
segInd{1} = ind;
load([cPath 'BehaviorVideo' filesep 'segInd2.mat'],'ind'); %load segment index
segInd{2} = ind;

%% recreate weights for video dimensions that were used in the linear model
if undoQR
    cBeta = NaN(size(redQRR,1),size(dimBeta,2));
    vidIdx = ismember(regIdx(~rejIdx), find(ismember(regLabels,{'bhvVideo'}))); % index for video regressors
    vidB = dimBeta(vidIdx, :);
    cInd = zeros(1,size(cBeta,1));
    cInd((end-sum(vidIdx))+1:end) = find(ismember(regLabels,{'bhvVideo'}));
    cBeta(cInd ~= 0,:) = vidB;

    % undo QR and equal length of video regressor columns
    betas = unQRVideo(cBeta, redQRR, cInd, regLabels)';
    betas = bsxfun(@rdivide,betas,sqrt(sum(betas.^2)));
else
    betas = dimBeta(regIdx(~rejIdx) == find(strcmp(regLabels, 'bhvVideo')), :)';
end

% un-zscore beta weights.
betas = bsxfun(@rdivide,betas,stdV(1:size(betas,2)));

% un-zscore vidV. This is later used to reconstruct tiled image as a sanity check.
vidV = bsxfun(@times,vidV,stdV);
vidV = bsxfun(@plus,vidV,meanV);
        
%% recreate and apply all temporal dimensions that were used in the segment SVDs    
%allocate reconstructed video
meanCam1 = NaN(numel(segInd{1}), 1, 'single');
meanCam2 = NaN(numel(segInd{2}), 1, 'single');
stdCam1 = NaN(numel(segInd{1}), 1, 'single');
stdCam2 = NaN(numel(segInd{2}), 1, 'single');
cam1 = NaN(numel(segInd{1}),size(betas,1), 'single');
orgCam1 = NaN(numel(segInd{1}),size(betas,1), 'single');
cam2 = NaN(numel(segInd{2}),size(betas,1),'single');
orgCam2 = NaN(numel(segInd{2}),size(betas,1),'single');

nBetas = betas * vidU(1:size(betas,2),:); %reconstructed betas for dimensions of first SVD. This is imagingDims x segmentVs
nV = vidV * vidU; %reconstruct temporal components from first SVD. This is time x segmentVs

segCnt = 0;
for iCams = 1:2
    for iSegs = 1:size(segInd{iCams},2)
    
        data = load([cPath 'BehaviorVideo' filesep 'SVD_Cam' int2str(iCams) '-Seg' int2str(iSegs) '.mat'],'U'); %load current segment
        cV = nV(:,segCnt + 1 : segCnt + size(data.U,1)); %data for current segment
        cV(isnan(cV(:,1)),:) = []; %remove NaNs

        if iCams == 1
            stdCam1(segInd{iCams}(:,iSegs),:) = sqrt(sum((data.U' * cov(cV)) .* data.U', 2));  %standard deviation map for current segment
            meanCam1(segInd{iCams}(:,iSegs),:) = (nanmean(cV,1) * data.U)';
            cam1(segInd{iCams}(:,iSegs),:) = (nBetas(:,segCnt + 1 : segCnt + size(data.U,1)) * data.U)';
        elseif iCams == 2
            stdCam2(segInd{iCams}(:,iSegs),:) = sqrt(sum((data.U' * cov(cV)) .* data.U', 2));  %standard deviation map for current segment
            meanCam2(segInd{iCams}(:,iSegs),:) = (nanmean(cV,1) * data.U)';
            cam2(segInd{iCams}(:,iSegs),:) = (nBetas(:,segCnt + 1 : segCnt + size(data.U,1)) * data.U)';            
        end
        segCnt = segCnt + size(data.U,1);
        
    end
end

%% reshape results to match original video resolution and save
movFiles = dir([cPath 'BehaviorVideo' filesep '*Video_*1.mj2']);
if ~isempty(movFiles)
    v = VideoReader([cPath 'BehaviorVideo' filesep movFiles(1).name]);
    temp = single(arrayResize(readFrame(v),2)); clear v
else
    temp = NaN(240, 320, 'single'); %if no video file is available, use standard resolution instead
end
meanCam1 = reshape(meanCam1,size(temp,1),size(temp,2),[]);
stdCam1 = reshape(stdCam1,size(temp,1),size(temp,2),[]);
cam1 = reshape(cam1,size(temp,1),size(temp,2),[]);

if ~isempty(movFiles)
    movFiles = dir([cPath 'BehaviorVideo' filesep '*Video_*2.mj2']);
    v = VideoReader([cPath 'BehaviorVideo' filesep movFiles(1).name]);
    temp = single(arrayResize(readFrame(v),2)); clear v
end
meanCam2 = reshape(meanCam2,size(temp,1),size(temp,2),[]);
stdCam2 = reshape(stdCam2,size(temp,1),size(temp,2),[]);
cam2 = reshape(cam2,size(temp,1),size(temp,2),[]);
clear temp

%save results
save([cPath fileExt 'betaBhvVid.mat'],'cam1','cam2','meanCam1','meanCam2','stdCam1','stdCam2');
