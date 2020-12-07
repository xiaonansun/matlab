
clear all;

rawVidDir = 'E:\uncompressed\FLEA3';

if rawVidDir(end) ~= filesep
    rawVidDir(end + 1) = filesep;
end

inVidExt = '.avi';
outVidExt = '.mj2';

dirContent = dir(fullfile([rawVidDir], ['*' inVidExt])); dirContent = {dirContent.name};

% fileIdx = 1;
% numOfFiles = 10;
numOfFiles = length(dirContent);

parfor fileIdx = 1:numOfFiles
    try
        inVidPath = [rawVidDir dirContent{fileIdx}];
        [fileDir, fileName, fileExt] = fileparts(inVidPath);
        outVidPath = fullfile(fileDir, [fileName outVidExt]);
        
        inVidObj = VideoReader(inVidPath);
        outVidObj = VideoWriter(outVidPath,'Archival');
        outVidObj.FrameRate = inVidObj.FrameRate;
        video = zeros(inVidObj.Height,inVidObj.Width,1,round(inVidObj.Duration*inVidObj.FrameRate));
        disp(['Converting file: ' [fileName fileExt]]);
        fIdx = 1;
        inVidObj.CurrentTime = 0;
        while hasFrame(inVidObj)
            video(:,:,1,fIdx) =rgb2gray(readFrame(inVidObj));
            fIdx = fIdx+1;
        end
        inVidObj.CurrentTime = 0;
        delete(inVidObj);
        
        open(outVidObj);
        writeVideo(outVidObj,uint8(video));
        close(outVidObj);
        disp(['Deleteing file: ' [fileName fileExt]]);
        delete(inVidPath);
    end

end

%% Loads mj2 into a matrix to compare with matrix imported from raw video file

mj2inVidPath = 'H:\mouse_single_stream_behavior_video_data\FLEA3\fc2_save_2019-04-22-162153-0000.mj2';
mj2VidObj = VideoReader(mj2inVidPath);
mj2video = zeros(mj2VidObj.Height,mj2VidObj.Width,1,round(mj2VidObj.Duration*mj2VidObj.FrameRate));
fIdx = 1;
mj2VidObj.CurrentTime = 0;
while hasFrame(mj2VidObj)
   mj2video(:,:,1,fIdx) =readFrame(mj2VidObj);
   fIdx = fIdx+1;
end

%%
diff = video-mj2video;
