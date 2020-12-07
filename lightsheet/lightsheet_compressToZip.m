function lightsheet_compressToZip(tiffDir)
% tiffDir = 'F:\WholeBrainImaging\AA1-PO16a';

% if ~exist('sPath','var') || isempty(sPath)
%     sPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
% end

if tiffDir(end) ~= filesep
    tiffDir(end + 1) = filesep;
end

inputExt = '*.tif';
outputExt = '.zip';

dirContent = dir(fullfile([tiffDir], fileExt)); dirContent = {dirContent.name};

parfor i = 1:length(dirContent)
    [tiffFileDir, tifFilename, tiffExt] = fileparts(dirContent{i});
    zipFilepath =[tiffDir filesep tifFilename outputExt];
    zip(zipFilepath,filenames);
    disp(['Compression file: ' dirContent{i}])
end

%%

