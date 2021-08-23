function Behavior_sampleVideoForLeap(cPath, nrClusters, nrSamples)

if ~exist('nrClusters','var') || isempty(nrClusters)
    nrClusters = 20;
end
if ~exist('nrSamples','var') || isempty(nrSamples)
    nrSamples = 500;
end
if cPath(end) ~= filesep
    cPath = [cPath filesep];
end
% cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\mSM63\SpatialDisc\01-Aug-2018\BehaviorVideo\';

sampleVideos = 25; %nr of video trials to sample from
overSample = 8; %nr of franes to use for SVD. This is multiplied with 'nrSamples'.
framesPerVideo = round(nrSamples * overSample / sampleVideos); %nr of frames to capture per video. will take from the end to have more licking frames.
nrDims = 50; %nr of dimensions for clustering

%% load combined dimensions and use kmeans to find different clusters
fprintf('Current recording: %s\nUsing %d clusters with %d samples each\n', cPath, nrClusters, nrSamples);
fprintf('Loading %d frames from %d video files\n', framesPerVideo, sampleVideos);

for iCams = 1 : 2
    tic
    %% get raw videos of current recording and perform clustering
    rawVids = dir([cPath '*_' num2str(iCams) '.mp4']);
    fileIdx = floor(length(rawVids)/sampleVideos) : floor(length(rawVids)/sampleVideos) : length(rawVids); %get frames from sample videos
    
    checker = true;
    for iFiles = fileIdx
        cFile = [cPath rawVids(iFiles).name];
        rawData = squeeze(importdata(cFile));
        rawData = squeeze(rawData(:,:,1,:));
        rawData = arrayResize(rawData, 2); %downsample by factor 2

        if checker
            allFrames = NaN(numel(rawData(:,:,1)), nrSamples * overSample);
            checker = false;
        end
        
        startIdx = find(isnan(allFrames(1,:)),1);
        
        if size(rawData, 3) > framesPerVideo
            rawData = rawData(:,:,end - framesPerVideo + 1 : end);
        end
        allFrames(:, startIdx : startIdx + size(rawData,3) - 1) = reshape(rawData,[], size(rawData,3));
    end
    allFrames(:, isnan(allFrames(1,:))) = []; %reject unused frames if not all of them were filled
    
    [~, s, vidV] = fsvd(allFrames, nrDims, 1, 0); %get svd
    vidV = s * vidV'; %multiply S into V, so only U and V from here on
    G = kmeans(vidV', nrClusters);
    clear vidV
    
    %% select random samples for each cluster and save
    idx = zeros(nrClusters, round(nrSamples / nrClusters));
    for iClust = 1 : nrClusters
        temp = find(G == iClust);
        if length(temp) > round(nrSamples / nrClusters)
            temp = temp(randperm(length(temp)));
            temp = temp(1:round(nrSamples / nrClusters));
        end
        idx(iClust, 1 : length(temp)) = temp;
    end
    idx = idx(:);
    idx(idx == 0) = [];
    
    allFrames = allFrames(:, idx);
    allFrames = reshape(allFrames, size(rawData,1), size(rawData,2), 1, size(allFrames,2));
    
    % write sample video for current cam
    if exist([cPath 'cam' num2str(iCams) '_sampleVideo.h5\'], 'file')
        delete([cPath 'cam' num2str(iCams) '_sampleVideo.h5\']);
    end
    h5create([cPath 'cam' num2str(iCams) '_sampleVideo.h5\'],'/box', size(allFrames),'ChunkSize', size(allFrames),'Datatype','uint8')
    h5write([cPath 'cam' num2str(iCams) '_sampleVideo.h5\'], '/box', uint8(allFrames))
    fprintf('Cam %d done. ', iCams); toc;
end
end
