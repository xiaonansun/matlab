function lightsheet_downscaleTiff(tiffDir)

% tiffDir = 'F:\WholeBrainImaging\AA1-PO16a';
% filePath = 'F:\WholeBrainImaging\AA1-PO16a\AA1-PO16a(0)-stitched_561_T001_Z001_C01.tif';

if tiffDir(end) ~= filesep
    tiffDir(end + 1) = filesep;
end
fileExt = '*.tif';
outputExt_j2k = '.j2k';
outputExt_jpg = '.jpg';
downFactor = 0.25;

dirContent = dir(fullfile([tiffDir], fileExt)); dirContent = {dirContent.name};

numOfFiles = length(dirContent);
% numOfFiles = 10;

parfor i = 1:numOfFiles
    [tiffFileDir, tifFilename, tiffExt] = fileparts(dirContent{i});
    tiffObj=Tiff([tiffDir filesep dirContent{i}],'r');
    imageData = read(tiffObj);
    imageData16d=imresize(imageData,downFactor,'bicubic');
    imwrite(imageData16d,[tiffDir filesep tifFilename '_downscaled_' num2str(1/downFactor) 'X' outputExt_j2k],'Mode','lossless');

    disp(['File ' num2str(i) '/' num2str(length(dirContent)) '. Downscaling file by ' num2str(1/downFactor) 'X: ' dirContent{i}]);
end


% imageData8=im2uint8(imageData);
% imwrite(imageData16d,[fileDir filesep fileName outputExt_jpg],'Quality',100);

% J = imresize(imageData,0.25);

% plotImg(imageData);