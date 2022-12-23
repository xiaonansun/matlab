% function lightsheet_convertJ2KToTiff(rootDir)

% This function reads converts JPEG2000 into Tiff files and then deletes the original Tiff files

% rootDir = 'D:\Churchland lab';

% fileDir = dir(rootDir); 
j2kDir = 'F:\LightsheetBatch2020-10-18\sample_14_tif';

tiffExt = '.tif';
j2kExt = '.j2k';

% numOfDir = length({fileDir.name});

% numOfFiles = 10;

% for idxDir = 1:numOfDir
%     tiffDir = [rootDir filesep fileDir(idxDir).name filesep]; 
dirContent = dir(fullfile(j2kDir, ['*' j2kExt])); 
dirContent = {dirContent.name};
    numOfFiles = length(dirContent);

    if isempty(dirContent)
        return
    else
        parfor idxJ2K = 1:numOfFiles
            try
            [j2kFileDir, j2kFilename, j2kExt] = fileparts(dirContent{idxJ2K});
            j2kPath = fullfile(j2kDir,[j2kFilename j2kExt]);
            tiffPath = fullfile(j2kDir,[j2kFilename tiffExt]);
%             tiffObj=Tiff(tiffPath,'r');
            img = imread(j2kPath);
            imwrite(img,tiffPath,'tif');
            disp(['Converting ' j2kPath ' to ' tiffPath]);
%             close(tiffObj);
%             delete(tiffPath);
            end
        end
    end

% end
