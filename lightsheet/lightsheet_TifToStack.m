function res = lightsheet_TifToStack(inputDir,outputDir,sampleID)

% This function invokes the custom scripts loadtiff.m and saveastiff.m

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

if isempty(dirContent)
   disp('No TIFF files exist in directory. Function terminated.');
   res = [];
   return
end

res = zeros(length(dirContent),1);

for i = 1:length(dirContent)
    tempImg = loadtiff(fullfile(dirContent(i).folder,dirContent(i).name));
    res(1) = saveastiff(tempImg,fullfile(outputDir,[sampleID '_5X.tif']),options);
    disp([num2str(i) 'TIFF(s) appended.']);
end