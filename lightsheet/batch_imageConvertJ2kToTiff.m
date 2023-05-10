function batch_imageConvertJ2kToTiff(inDir)

%% 

% inDir = '/Users/xiaonansun/Documents/stp_data/170622_KM_Fezfh2bGFPmale_processed/stitchedImage_ch1'; % comment this line if not prototyping

% This function converts lossy JPEG2000 (.jp2) files into tiff files while
% preserving the bit depth.

% XnViewMP unfortunately cannot perform this conversion on STP data from
% the repository (brainimagelibrary.org) without distorting and/or reducing
% the bit depth.
% IrfanView 32-bit can perform this function but cannot be used on MacOS

% To send this command as a batch for processing use:
% j = batch(@imageConvertJ2kToTiff,0,{'inDir'})

outDir = inDir;

outputFileExt = '.tif';
j2kExt = '.jp2';

dirContent = dir(fullfile(inDir, ['*' j2kExt])); dirContent = {dirContent.name};

tic
parfor idx = 1:length(dirContent)
        disp(['Working on file: ' dirContent{idx}]);
        [~, imgFilename, ~] = fileparts(dirContent{idx});
        imgPath = fullfile(inDir,dirContent{idx});
        imgMat = imread(imgPath);
        imwrite(imgMat, fullfile(outDir, [imgFilename outputFileExt]));
end

fprintf(['Time completed: ' char(datetime) '.\nConversion completed in ' num2str(toc) ' seconds.\nAll files saved in ' outDir '\n']);

if ismac
disp('Open folder <a href="matlab: finder(outDir), ">here</a>.'); % this line only works on MacOS
end


