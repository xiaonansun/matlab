function data = twoP_zscoreTrializeContinuousSDU(animal, session, overwriteSDU)

%% Trializes continuous zscore traces that are already computed from 
% animal = 'Fez61'; session = '20200604'; overwriteSDU = true;

S = twoP_settings;

if ~exist('overwriteSDU','var') || isempty(overwriteSDU)
    overwriteSDU = true;
end

% Load neural data struct (data.mat)
neural_filepath = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'data.mat');
data = load(neural_filepath); data = data.data;
if isfield(data,'sdu') && overwriteSDU == 1
    disp('data.sdu exists. Will overwrite.')
elseif isfield(data,'sdu') && overwriteSDU == 0
    disp('data.sdu exists. Will not overwrite.')
    sdu = data.sdu;
    return
end

% Load SDU data
sdu_files = dir(fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'zscore*.mat'));
if length(sdu_files) > 1
    error('There are multiple zscored traces, please specify which one to process.');
elseif isempty(sdu_files)
    error('There are no z-scored traces, please first compute zscored traces from spks.');
elseif length(sdu_files) == 1
    sdu = load(fullfile(sdu_files.folder,sdu_files.name)); sdu = sdu.sdu;
end

isneuron_filepath = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'iscell.npy');
isneuron = readNPY(isneuron_filepath); isneuron = logical(isneuron(:,1));

%%
% Define stimulus-onset frames
if iscell(data.stimFramesOrig)
    stimFrames = data.stimFramesOrig;
    warning('data.stimFramesOrig is a cell, imaging was interrupted, will combine cells into matrix.')
else
    stimFrames = data.stimFramesOrig;
end

%% Loop through all frames including data structs organized as cells
if iscell(stimFrames)
    for iC = 1:length(stimFrames)
        if length(stimFrames{iC}) ~= size(data.neural{iC},3)
            error('The number stimulus-on frame indices does not equal the number recorded frames: check length(data.stimFramesOrig) against size(data.neural,3)');
        end
        data.sdu{iC} = zeros(size(data.neural{iC},1),size(data.neural{iC},2),size(data.neural{iC},3));
        nTrials = length(stimFrames{iC});
        for tr = 1:nTrials
            try
                % PROBLEM: stimFrames is the frame count of 
                data.sdu{iC}(:, :, tr) = sdu{iC}(isneuron, stimFrames{iC}(tr) + data.neuralTimes{iC});
            catch ME
                disp(ME.message);
            end
        end
    end
else
    if length(stimFrames) ~= size(data.neural,3)
        error('The number stimulus-on frame indices does not equal the number recorded frames: check length(data.stimFramesOrig) against size(data.neural,3)');
    end
    data.sdu = zeros(size(data.neural,1),size(data.neural,2),size(data.neural,3));
    nTrials = length(stimFrames);
    for tr = 1:nTrials
        try
            data.sdu(:, :, tr) = sdu(isneuron, stimFrames(tr) + data.neuralTimes);
        catch ME
            disp(ME.message);
        end
    end
end
% sdu = data.sdu;

save(neural_filepath,'data');

disp(['Z-scored traces of ' animal ' ' session ' is trialized and saved as ' neural_filepath])
