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

if ~exist('sampleID','var') || isempty(sampleID)
    sampleID = regexp(inputDir,'\w*sample\w*','match');
    if isempty(sampleID)
    disp('Directory path of tiff files missing: user must specify input directory (inputDir)! Function terminated.');
    return
    end
end

options.big = true;
options.append = true;
options.message = false;

dirContent = dir(fullfile(inputDir,'*X.tif'));
dirContent = natsortfiles({dirContent.name});
idxFNtail = strfind(dirContent{1},'X.tif');
saveFileName = fullfile(outputDir,[sampleID{:} '_' dirContent{1}(idxFNtail-1) 'X.tif']);

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