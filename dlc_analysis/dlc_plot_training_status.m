function dlc_plot_training_status(projectDir)

%  projectDir is the root directory of the DLC project(s) 
pDir=projectDir;
csvfName='learning_stats.csv';
pDirList = split(genpath(pDir),';');
csvFileList = cell(1,length(pDirList));

for i = 1:length(pDirList)
    csvFileList(i) = dir(fullfile(pDirList{i},csvFName));
end