clear all;
vidDir = 'Z:\BehaviorVideo\mSM85\SpatialDisc\Session Data\mSM85_SpatialDisc_Jan22_2020_Session1';
S = regexp(vidDir,filesep,'split');
face = '_1'; mouth = '_2'; ext = '.mp4';

faceFiles = dir(fullfile(vidDir,['*' face ext]));
mouthFiles = dir(fullfile(vidDir,['*' mouth ext]));

faceList = {faceFiles.name}';   
mouthList = {mouthFiles.name}';   

faceOutputVideo = VideoWriter(fullfile(vidDir,[S{end} face]),'Motion JPEG 2000'); 
faceOutputVideo.CompressionRatio = 10;
% faceOutputVideo.Quality = 75;
mouthOutputVideo = VideoWriter(fullfile(vidDir,[S{end} mouth]),'Motion JPEG 2000'); 
mouthOutputVideo.CompressionRatio = 10;
% mouthOutputVideo.Quality = 75;
% if all clips are from the same source/have the same specifications
% just initialize with the settings of the first video in videoList
faceInputVideo_init = VideoReader([vidDir filesep faceList{1}]); % first video
% faceInputVideoMat = zeros(faceInputVideo_init.Height, faceInputVideo_init.Width,1,faceInputVideo_init.Duration*faceInputVideo_init.FrameRate);
mouthInputVideo_init = VideoReader([vidDir filesep mouthList{1}]); % first video

faceOutputVideo.FrameRate = faceInputVideo_init.FrameRate;
mouthOutputVideo.FrameRate = mouthInputVideo_init.FrameRate;

open(faceOutputVideo); open(mouthOutputVideo);

% numFaceVids = 10;
% numMouthVids = 10;
numFaceVids = length(faceList);
numMouthVids = length(mouthList);

for i = 1:numFaceVids
    faceInputVideo = VideoReader([vidDir filesep faceList{i}]);
    mouthInputVideo = VideoReader([vidDir filesep mouthList{i}]);
    faceInputVideoMat = zeros(faceInputVideo.Height, faceInputVideo.Width,1,ceil(faceInputVideo.Duration*faceInputVideo.FrameRate));
    mouthInputVideoMat = zeros(mouthInputVideo.Height, mouthInputVideo.Width,1,ceil(mouthInputVideo.Duration*mouthInputVideo.FrameRate));
    fIdx = 1;
    mIdx = 1;
    while hasFrame(faceInputVideo)
        faceInputVideoMat(:,:,1,fIdx) = rgb2gray(readFrame(faceInputVideo));
        fIdx = fIdx+1;
    end
    while hasFrame(mouthInputVideo)
        mouthInputVideoMat(:,:,1,mIdx) = rgb2gray(readFrame(mouthInputVideo));
        mIdx = mIdx+1;
    end
    disp(faceList{i});
    disp(mouthList{i});
    
    writeVideo(faceOutputVideo, uint16(faceInputVideoMat));
    writeVideo(mouthOutputVideo, uint16(mouthInputVideoMat));
end
close(faceOutputVideo);
close(mouthOutputVideo);
