function cmbBhv = twoP_RepairSession(animal,cell_of_filenames,trial_codes)

%% combine behavior
% basePath = 'G:\Google Drive\Behavior_Simon\PLEX51\SpatialDisc\Session Data\';
S = twoP_settings;

% basePath = 'Y:\data\Behavior_Simon\Plex51\SpatialDisc\Session Data';
if length(cell_of_filenames) < 2
    disp('There are less than two behavior files for this session. No need to append behavior files.');
    return
elseif length(cell_of_filenames) == 2
    disp('There are two behavior files for this session, will append these files.')
    firstPath = fullfile(S.dir.bhvRootDir,animal,S.dir.bhvSubDir,[cell_of_filenames{1} '.mat']);
    secondPath = fullfile(S.dir.bhvRootDir,animal,S.dir.bhvSubDir,[cell_of_filenames{2} '.mat']);
end

% firstPath = [basePath 'Plex51_SpatialDisc_Apr01_2020_Session1'];
% secondPath = [basePath 'Plex51_SpatialDisc_Apr01_2020_Session2'];

% gapTrials = 1; %nr of imaging trials that were recorded between behavioral sessions

% [a, b] = fileparts(firstPath);
% load([a filesep b '.mat'])
load(firstPath,'SessionData');
idx1 = 1:length(SessionData.Rewarded); %change this if only a subset of data should be used
bhv1 = selectBehaviorTrials(SessionData,idx1);

idxPk = findpeaks(trial_codes);
gapTrials = trial_codes(idxPk)-length(bhv1.Rewarded); %nr of imaging trials that were recorded between behavioral sessions

gapBhv = selectBehaviorTrials(SessionData, 1:gapTrials); %add empty trials (this is to compensate trials where bpod broke)
gapBhv.DidNotLever = true;
gapBhv.DidNotChoose = true;


bhv1g = appendBehavior(bhv1,gapBhv); %combine bhv data

% [a, b] = fileparts(secondPath);
% load([a filesep b '.mat'])
load(secondPath,'SessionData');
idx2 = 1:length(SessionData.Rewarded); %change this if only a subset of data should be used
bhv2 = selectBehaviorTrials(SessionData,idx2);

cmbBhv = appendBehavior(bhv1g,bhv2); %combine bhv data
cmbBhv.nTrials = idx1(end)+idx2(end)+gapTrials;
cmbBhv.bhv1 = bhv1; cmbBhv.bhv2 = bhv2; cmbBhv.gapBhv = gapBhv;
% temp = findstr(b,'_');
% cFile = [a filesep b(1:temp(end)) 'Combined'];
% save(cFile, 'SessionData');

%% Rename behavior video to get matching file names - this one is needed if you want to use behavioral videos.
% Widefield_RenameBehaviorVideo(firstPath, secondPath, idx1(end)+gapTrials)
