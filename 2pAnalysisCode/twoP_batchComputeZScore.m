S = twoP_settings;

findFileName = 'spks.npy';

filelist = dir(fullfile(S.dir.imagingRootDir,['**\' findFileName]));
folderList = {filelist.folder};
splitCells = cellfun(@(x) regexp(x,filesep,'split'),folderList,'UniformOutput',false);
animalID = cellfun(@(x) x{7},splitCells,'UniformOutput',false); % Parsed animalID from the directory list
sessionID = cellfun(@(x) x{9},splitCells,'UniformOutput',false); % parsed sessionID from the directory list

%%
parfor i = 1:length(animalID)
    twoP_computeZScoreFromSpks(animalID{i},sessionID{i},1);
end

%%
S = twoP_settings;

findFileName = 'data.mat';

filelist = dir(fullfile(S.dir.imagingRootDir,['**\' findFileName]));
folderList = {filelist.folder};
splitCells = cellfun(@(x) regexp(x,filesep,'split'),folderList,'UniformOutput',false);
animalID = cellfun(@(x) x{7},splitCells,'UniformOutput',false); % Parsed animalID from the directory list
sessionID = cellfun(@(x) x{9},splitCells,'UniformOutput',false); % parsed sessionID from the directory list


for i = 1:length(animalID)
    twoP_alignZScore(animalID{i},sessionID{i});
end
