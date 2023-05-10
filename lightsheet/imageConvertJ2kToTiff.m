function res = imageConvertJ2kToTiff(sampleID)
% This function opens J2K files, re-scales the images, and then saves as
% tiff files

% To send this command as a batch for processing use:
% j = batch(@lightsheet_downscale_convert_J2kToTiff,0,{'sample_N_tif'})

baseDir = 'F:\LightsheetBatch2020-10-18';
imgDir = fullfile(baseDir,sampleID);
dsImgDir = fullfile(imgDir,'downscaled');
dsImgDirCont = dir(dsImgDir); dsImgDirCont = {dsImgDirCont.name};
fileExt = '*.tif';
j2kExt = '.j2k';

dFactor = 0.2;

dirContent = dir(fullfile(imgDir, ['*' j2kExt])); dirContent = {dirContent.name};

if ~exist(dsImgDir,'dir') || isempty(dsImgDir) % creates a new directory with the downscaled files
    mkdir(dsImgDir);
end

tic
parfor idx = 1:length(dirContent)
    try
        disp(['Working on file: ' dirContent{idx}]);
        [~, imgFilename, ~] = fileparts(dirContent{idx});
        imgPath = fullfile(imgDir,dirContent{idx});
        imgMat = imresize(imread(imgPath),0.25,'bicubic');
        imwrite(imgMat, fullfile(dsImgDir, [imgFilename 'downscaled_' num2str(round(1/dFactor)) 'X.tif']));
    end
end
disp(['Conversion completed in ' num2str(toc) 'seconds.']);


