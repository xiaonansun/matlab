function res = lightsheet_downscale_convert_J2kToTiff(sampleID)


% This function reads converts Tiff into JPEG2000 files and then deletes the original Tiff files
% Conversion is lossless

% To send this command as a batch for processing use:
% j = batch(@lightsheet_downscale_convert_J2kToTiff,0,{'sample_N_tif'})

% rootDir = 'G:';
% subDir = dir(rootDir);
% subDir = subDir(~ismember({subDir.name},{'.','..'}));
% sampleID = 'sample_10_tif';
baseDir = 'F:\LightsheetBatch2020-10-18';
imgDir = fullfile(baseDir,sampleID);
dsImgDir = fullfile(imgDir,'downscaled');
dsImgDirCont = dir(dsImgDir); dsImgDirCont = {dsImgDirCont.name};
fileExt = '*.tif';
j2kExt = '.j2k';
mj2Ext = '.mj2';

dFactor = 0.2;

dirContent = dir(fullfile(imgDir, ['*' j2kExt])); dirContent = {dirContent.name};

if ~exist(dsImgDir,'dir') || isempty(dsImgDir) % creates a new directory with the downscaled files
    mkdir(dsImgDir);
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
    res = lightsheet_TifToStack(dsImgDir,[],sampleID);
    disp('Downscaled files are combined into a single tiff');
else
    dsDirCont = dir(fullfile(dsImgDir,'*downscaled_5X.tif'));
    dsDirCont = {dsDirCont.name};
    combFilePath = fullfile(dsImgDir,[sampleID '_5X.tif']);
    combFileExist = ~isempty(dir(combFilePath));
    res = lightsheet_TifToStack(dsImgDir,[],sampleID);
    disp('Downscaled files are combined into a single tiff');
end



% numOfFiles = length(dirContent);

%% SAVE AS STACK
% use the loadtiff and saveastiff functions
% imgDir = fullfile(baseDir,sampleID); S = regexp(imgDir,filesep,'split');
%
% tiffDir = fullfile(imgDir,'downscaled');
%
% res = lightsheet_TifToStack(dsImgDir,[],sampleID);