
% function lightsheet_downscale_convert_J2kToTiff(tiffDir)

% This function reads converts Tiff into JPEG2000 files and then deletes the original Tiff files
% Conversion is lossless

clear all;


% rootDir = 'G:';
% subDir = dir(rootDir); 
% subDir = subDir(~ismember({subDir.name},{'.','..'}));
imgDir = 'W:\lightsheet\batch03\sample_2'; S = regexp(imgDir,filesep,'split');

if imgDir(end) ~= filesep
    imgDir(end + 1) = filesep;
end

if ~exist([imgDir 'downscaled']) % creates a new directory with the downscaled files
    mkdir([imgDir 'downscaled']);
end

fileExt = '*.tif';
j2kExt = '.j2k';
mj2Ext = '.mj2';

dFactor = 0.25;

dirContent = dir(fullfile(imgDir, ['*' j2kExt])); dirContent = {dirContent.name};

numOfFiles = length(dirContent);
% numOfFiles = 10;  % for testing purposes

tic
parfor idx = 1:numOfFiles
    try
        disp(['Working on file: ' dirContent{idx}]);
        [imgFileDir, imgFilename, imgExt] = fileparts(dirContent{idx});
        imgPath = [imgDir dirContent{idx}];
        imgMat = imresize(imread(imgPath),0.25,'bicubic');
        imwrite(imgMat, fullfile(imgDir, 'downscaled', [imgFilename 'downscaled_' num2str(round(1/dFactor)) 'X.tif']));
%         writeVideo(vObj, uint16(imgPath));
    end
end
disp(['Conversion completed in ' num2str(toc) 'seconds.']);

