function res = lightsheet_TifToStack(inputDir,outputDir,sampleID)

% This function invokes the custom scripts loadtiff.m and saveastiff.m
% To send this function as a batch, use: 
% j1 = batch(@lightsheet_TifToStack,1,{'F:\LightsheetBatch2020-10-18\sample_N_tif\downscaled',[],'sample_N_tif'})

if ~exist('inputDir','var') || isempty(inputDir)
    disp('Directory path of tiff files missing: user must specify input directory (inputDir)! Function terminated.');
    return
end

if ~exist('outputDir','var') || isempty(outputDir)
    outputDir = inputDir;
end

options.big = true;
options.append = true;
options.message = false;

dirContent = dir(fullfile(inputDir,'*5X.tif'));
dirContent = natsortfiles({dirContent.name});
saveFileName = fullfile(outputDir,[sampleID '_5X.tif']);

if exist(saveFileName,'file')
    delete(saveFileName);
end

if isempty(dirContent)
   disp('No TIFF files exist in directory. Function terminated.');
   res = [];
   return
end

res = zeros(length(dirContent),1);

for i = 1:length(dirContent)
    tempImg = loadtiff(fullfile(inputDir,dirContent{i}));
    res(1) = saveastiff(tempImg,saveFileName,options);
    disp([num2str(i) 'TIFF(s) appended.']);
end