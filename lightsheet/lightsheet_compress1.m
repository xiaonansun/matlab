tiffDir = 'E:\2P_data_shared\test_data\compress';
if tiffDir(end) ~= filesep
    tiffDir(end + 1) = filesep;
end

subDir = [tiffDir 'img_test'] ;
if ~exist(subDir, 'dir')
    mkdir(subDir);
end

fileExt = '*.tif';
vidExt = '.mj2';
outputExt_j2k = '.j2k';

dirContent = dir(fullfile([tiffDir], fileExt)); dirContent = {dirContent.name};

idxTiff = 1;

tiffPath = [tiffDir dirContent{idxTiff}];
[tiffFileDir, tifFilename, tiffExt] = fileparts(tiffPath);

imgStack = [];
tifInfo = imfinfo(tiffPath);
numOfImages = length(tifInfo); 

tic
parfor k = 1:numOfImages
    tifImg=imread(tiffPath, k);
    imgStack(:,:,1,k)=tifImg;
    imwrite(tifImg,[subDir filesep tifFilename '_' num2str(k) outputExt_j2k],'Mode','lossless');
end
toc

zip([tiffDir filesep tifFilename '.zip'],[tiffDir filesep tifFilename fileExt(2:end)]);

vidObj = VideoWriter([tiffFileDir filesep tifFilename vidExt],'Archival');
open(vidObj);
writeVideo(vidObj,imgStack);
% pngPath = [tiffDir tifFilename outputExt_png];

tiffObj=Tiff(tiffPath,'r');
tiff = read(tiffObj);