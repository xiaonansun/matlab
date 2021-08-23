function [npy,data,bhv,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(varargin)

% Inputs:
% (1) animal: the name of the animal, e.g. 'CSP27'
% (2) session: name of the session, e.g. '200329'
% (3) doSmooth: boolean
% (4) filterWindow: savitsky golay filter windows length
% (5) 

% Outputs
% (1) npy: 
% (2) data: two photon data in a 3D matrix
% (3) bhv: behavior data (SessionData)
% (4) bhvFilePath: path of the Bpod behavior file
% (5) suite2pDir: directory of the s2p-processed 2p data files


animal = varargin{1};
session = varargin{2};

if nargin > 2
    doSmooth = varargin{3};
    filterWindow = varargin{4};
    loadSavedData = varargin{5};
else
    doSmooth = true;
    filterWindow = 10;
    loadSavedData = false;
end

analysisFileName = [animal '_' session];

% Define directories

D = struct;
if convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('MANHASSET')
    imagingRootDir = 'G:\2PData';
elseif convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('SUNHP')
%     imagingRootDir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
    imagingRootDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
    codeRootDir = 'C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab';
    disp('Loading suite2p data...');
elseif convertCharsToStrings(computer) == convertCharsToStrings('GLNXA64') || isunix == 1
    imagingRootDir = '/grid/churchland/data/data/richard_s2p_npy';
    codeRootDir = '/grid/churchland/home/xisun/matlab';
end

D.imagingRootDir = imagingRootDir; D.codeRootDir = codeRootDir;
addpath(genpath(codeRootDir)); % adds the paths of the code files to be executed

imagingSubDir = ['suite2p' filesep 'plane0'];

disp(['Animal: ' animal '; Session: ' session]);

% % --- Load the google sheets document "2photon acquisition record" --- %
% docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';
% expTable=GetGoogleSpreadsheet(docid); % this function (GetGoogleSpreadsheet.m) needs to be downloaded
% bhvColIdx=find(contains(expTable(1,:),'Behavior file name'));
% iFolderColIdx=find(contains(expTable(1,:),'Folder'));
% 
% try
% bhvRowIdx = find(contains(expTable(:,iFolderColIdx),session));
% bhvFName = expTable{bhvRowIdx(contains(expTable(bhvRowIdx),animal)),bhvColIdx};
% catch ME
%     disp([ME.message]);
%     disp('Cannot load session. Please check the session name input.');
%     analysis.error.behaviorTable = ME.message;
% end

analysisDir = fullfile(imagingRootDir,'analysis'); %directory of the analysis history summary file

if ~exist(fullfile(imagingRootDir,'analysis'),'dir')
    mkdir(analysisDir);
end

try
load(fullfile(analysisDir,'analysis.mat'),'analysis'); % load the analysis history summary file
currentDate=datestr(now,'yyyy-mm-dd_HHMMss'); currAnalIdx=length(analysis)+1; % append data of the current analysis to the summary file

fprintf('Loading 2P imaging data...');

% organize suite2p data into npy files (ROI and spiking)
tic

s2pDataDir = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir);
codeDir = fullfile(codeRootDir,'2pAnalysisCode');
npy = twoP_importSuite2pData(animal,session, s2pDataDir, codeDir); 
disp(['The .npy file was loaded in ' num2str(toc) ' seconds.'])

% IMPORTANT!!! smoothes inferred spikes with a gaussian filter
if doSmooth==true
    spks = smoothCol(npy.spks,2,filterWindow,'gauss');
else
    spks = npy.spks;
end

tic
try
data = twoP_alignDetectionTask(npy.ops, spks, npy.iscell, npy.redcell, npy.bin_MScan_filepath); % align suite2p data to sensory stimulus
catch ME
end
disp(['Trial alignment to stimulus onset was completed in ' num2str(toc) ' seconds.'])

if exist('ME','var')
    disp(['Error detected: ' ME.identifier '. ' ME.message])
else
    disp('2P DATA LOADED!');
end

% Load behavior data
fprintf('Loading behavior data...');
[bhv,bhvFilePath] = twoP_loadBehaviorSession(animal,session); 
% bhv = bhv.bhv;
% [filepath,bhvFileName,ext]=fileparts(bhvFilePath); clear filepath ext;
fprintf('DONE!\n');

suite2pDir = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir); data.suite2pDir = suite2pDir;
disp(['Directory of suite2p output: ' suite2pDir]);
disp(['Path of behavior data: ' bhvFilePath]);
disp(['Number of trial codes received by MScan: ' num2str(max(data.trialNumbers))]);
disp(['Number of trials included in neural data matrix (MScan analog data not rejected): ' num2str(length(data.trialNumbers))]);
disp(['Number of trials recorded by Bpod: ' num2str(length(bhv.Rewarded))]);
if abs(max(data.trialNumbers) -length(bhv.Rewarded)) >= 10; disp('Behavior and imaging differs by more than 10 trials: do you have the correct behavior file?'); end

analysis(currAnalIdx).date=currentDate;
analysis(currAnalIdx).animal=animal; analysis(currAnalIdx).session=session;
if exist('ME','var')
analysis(currAnalIdx).errors=ME.message;
end
analysis(currAnalIdx).twopRecordedTrials=max(data.trialNumbers);
analysis(currAnalIdx).twopAlignedDataTrials = length(data.trialNumbers);
analysis(currAnalIdx).bpodTrials = length(bhv.Rewarded);
save([analysisDir filesep 'analysis.mat'],'analysis'); % save a historical record of this analysis script that has been run
writetable(struct2table(analysis),[analysisDir filesep 'analysis.csv'])

end