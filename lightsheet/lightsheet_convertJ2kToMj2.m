% function lightsheet_compressTiffToJ2K(tiffDir)

% This function reads converts Tiff into JPEG2000 files and then deletes the original Tiff files
% Conversion is lossless

clear all;


% rootDir = 'G:';
% subDir = dir(rootDir); 
% subDir = subDir(~ismember({subDir.name},{'.','..'}));
imgDir = 'F:\LightsheetBatch2020-10-18\sample_14_tif'; S = regexp(imgDir,filesep,'split');
vFileName = S{end};

if imgDir(end) ~= filesep
    imgDir(end + 1) = filesep;
end

fileExt = '*.tif';
j2kExt = '.j2k';
mj2Ext = '.mj2';

dFactor = 0.25;

% imgDir = [rootDir filesep subDir(idxDir).name filesep];
dirContent = dir(fullfile(imgDir, ['*' j2kExt])); dirContent = {dirContent.name};

numOfFiles = length(dirContent);
% numOfFiles = 2;

% for idxDir = 1:numOfDir
%     tiffDir = [rootDir filesep subDir(idxDir).name filesep]; 
%     dirContent = dir(fullfile(tiffDir, fileExt)); dirContent = {dirContent.name};
%     numOfFiles = length(dirContent);
% 
%     if isempty(dirContent)
%         continue
%     else
vObj = VideoWriter([imgDir vFileName],'Motion JPEG 2000');
vObj.FrameRate = 10;
open(vObj);

img1=imread([imgDir dirContent{1}]);
img = zeros(round(size(img1,1)*dFactor),round(size(img1,2)*dFactor),1,numOfFiles);
clear img1;

for idx = 1:numOfFiles
    try
        disp(['Working on file: ' dirContent{idx}]);
        [imgFileDir, imgFilename, imgExt] = fileparts(dirContent{idx});
        imgPath = [imgDir dirContent{idx}];
        currentImage = imresize(imread(imgPath),0.25,'bicubic');
        writeVideo(vObj, uint16(imgPath));
    end
end

% writeVideo(vObj,img);
close(vObj);
%     end
% 
% end
