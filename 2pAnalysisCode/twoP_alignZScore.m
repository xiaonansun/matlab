function twoP_alignZScore(animal,session)
%%
% animal = 'CSP27'; session = '20200319';
file_name_prefix = 'zscore';
S = twoP_settings;
load_dir = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir);
load_file = dir(fullfile(load_dir,[file_name_prefix '*']));
load_file = {load_file.name};

if length(load_file) > 1
    disp('There are more than one file that can be loaded, please change variable file_name_prefix');
    return
end

load(fullfile(load_dir,load_file{:}));

if exist('zscore','var')
    disp(['The z-score trace has already been aligned for ' animal ' ' session '.']);
    return
end

load(fullfile(load_dir,'data.mat'));

stimFrames = data.stimFramesOrig; % frame index of stimulus onset
if iscell(stimFrames)
    load(fullfile(load_dir,'ops.mat')); ops = ops{1};
    idx_reset_file = find(diff(ops.frames_per_file) > 1);
    for i = 1:length(idx_reset_file)
        idx_reset(i) = sum(ops.frames_per_file(1:idx_reset_file(i)-1))+1;
        stimFrames{i+1}=stimFrames{i+1}+sum(ops.frames_per_file(1:idx_reset_file(i)-1));
    end
    tmp = [];
    for i = 1:length(stimFrames)
        tmp = [tmp stimFrames{i}];
    end
    stimFrames = tmp;
end

idxRelStim = data.neuralTimes; % epoch frame indices relative to stim onset
if iscell(idxRelStim)
    idxRelStim = idxRelStim{1};
end

isneuron = readNPY(fullfile(load_dir,'iscell.npy'));

sdu_iscell = sdu(logical(isneuron(:,1)),:);
zscore = zeros(sum(isneuron(:,1)),length(idxRelStim),length(stimFrames)); 
for i = 1:length(stimFrames)
    zscore(:,:,i) = sdu_iscell(:,stimFrames(i) + idxRelStim);
end
zscore = single(zscore); 

save(fullfile(load_dir,load_file{:}),'zscore','sdu','-v7.3');
disp(['The z-score trace has been aligned for ' animal ' ' session '.']);
%%
iCell = 5; tr = 60;
plot(zscore(iCell,:,tr),data.neural(iCell,:,tr),'.k');