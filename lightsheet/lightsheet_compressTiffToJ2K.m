function lightsheet_compressTiffToJ2K(rootDir)

% This function reads converts Tiff into JPEG2000 files and then deletes the original Tiff files
% Conversion is lossless

% rootDir = 'D:\Churchland lab';
subDir = dir(rootDir); 
subDir = subDir(~ismember({subDir.name},{'.','..'}));

% if tiffDir(end) ~= filesep
%     tiffDir(end + 1) = filesep;
% end

fileExt = '*.tif';
j2kExt = '.j2k';

numOfDir = length({subDir.name});

% numOfFiles = 10;

for idxDir = 1:numOfDir
    tiffDir = [rootDir filesep subDir(idxDir).name filesep]; 
    dirContent = dir(fullfile(tiffDir, fileExt)); dirContent = {dirContent.name};
    numOfFiles = length(dirContent);

    if isempty(dirContent)
        continue
    else
        parfor idxTiff = 1:numOfFiles
            try
            [tiffFileDir, tifFilename, tiffExt] = fileparts(dirContent{idxTiff});
            tiffPath = [tiffDir dirContent{idxTiff}];
            jpgPath = [tiffDir tifFilename j2kExt];
            tiffObj=Tiff(tiffPath,'r');
            tiff = read(tiffObj);
            imwrite(tiff,jpgPath,'Mode','lossless');
            disp(['Working on file: ' dirContent{idxTiff}]);
            close(tiffObj);
            delete(tiffPath);
            end
        end
    end

end
