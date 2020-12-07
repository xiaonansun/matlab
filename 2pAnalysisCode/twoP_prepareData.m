function twoP_prepareData(tPath,sPath,trialIdx)

if ~strcmpi(tPath(end),filesep)
    tPath = [tPath filesep];
end
if ~strcmpi(sPath(end),filesep)
    sPath = [sPath filesep];
end

%% get 2p data and trial index from behavior
datFile = dir([sPath '*.mat']);
binFile = dir([sPath '*.bin']);
data = twoP_alignDetectionTask([sPath datFile.name], [sPath binFile.name], true);
twoP_checkTrialDrift(data, 0.4, 0.6); %check for drift over trials

bTrials = twoP_FindGoodTrials(tPath);
if ~isinf(trialIdx)
    bTrials(~ismember(bTrials,data.trialNumbers(trialIdx))) = []; %reject trials based on pre-determiend trial index. This is in case something bad happened in the session.
end
data.bhvTrials = bTrials;

% save 2p data to new folder
save([tPath 'data.mat'], 'data', 'bTrials');

% move textfiles
txtFiles = dir([sPath '*.txt']);
for iFiles = 1:length(txtFiles)
    copyfile([sPath txtFiles(iFiles).name], [tPath txtFiles(iFiles).name]);
end
    
