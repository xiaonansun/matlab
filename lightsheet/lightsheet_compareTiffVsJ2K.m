tiffDir = 'F:\WholeBrainImaging\AA1-PO16a';

if tiffDir(end) ~= filesep
    tiffDir(end + 1) = filesep;
end
fileExt = '*.tif';
outputExt_j2k = '.j2k';

dirContent = dir(fullfile([tiffDir], fileExt)); dirContent = {dirContent.name};

% idxTiff = 1500;
% numOfFiles = length(dirContent);
numOfFiles = 1;

maxmin = zeros(numOfFiles,2);

parfor idxTiff = 1:numOfFiles
[tiffFileDir, tifFilename, tiffExt] = fileparts(dirContent{idxTiff});
tiffPath = [tiffDir dirContent{idxTiff}];
pngPath = [tiffDir tifFilename outputExt_j2k];
tiffObj=Tiff(tiffPath,'r');
tiff = read(tiffObj);
imwrite(tiff,pngPath,'Mode','lossless');
png = imread(pngPath);
delta = tiff-png;
maxmin(idxTiff,:) = [max(delta(:)) min(delta(:))];
disp(['Working on file: ' dirContent{idxTiff}]);
end

% histogram(tiff);
% histogram(jpg);
