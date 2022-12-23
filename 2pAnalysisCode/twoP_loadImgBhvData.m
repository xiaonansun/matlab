function [data,bhv,npy]=twoP_loadImgBhvData(varargin)

%% Loads suite2p-processed two-photon imaging data and behavior data
% and trializes imaging data by aligning to the onset of the 
% sensory stimulus

% Inputs:
% (1) animal: the name of the animal, e.g. 'CSP27'
% (2) session: name of the session, e.g. '200329'
% (3) doSmooth: boolean (1 for yes, 0 for no)
% (4) filterWindow: filter windows length
% (5) loadSavedData: load data that is already saved as data.mat

% Outputs
% (1) npy:
% (2) data: two photon data in a 3D matrix (neuron x frames x trials)
% (3) bhv: behavior data (SessionData)

disp(newline);

% First load some constants and required data
exps = twoP_getAcquisitionRecord; % loads the 2p acquisition record from the google sheets file
S = twoP_settings;

% Specify the animal and the session IDs
animal = varargin{1};
session = varargin{2};

if nargin > 2
    doSmooth = varargin{3};
    if isempty(doSmooth); doSmooth = false; end
    filterWindow = varargin{4};
    if isempty(filterWindow); filterWindow = 10; end
    loadSavedData = varargin{5};
    if isempty(loadSavedData); loadSavedData=false; end
else
    doSmooth = true;
    filterWindow = 10;
    loadSavedData = false;
end

% analysisFileName = [animal '_' session];

% Define working directories 
imagingRootDir = S.dir.imagingRootDir;
codeRootDir = S.dir.codeRootDir;
bhvRootDir = S.dir.bhvRootDir;
imagingSubDir = S.dir.imagingSubDir;
bhvSubDir = S.dir.bhvSubDir;
s2pDataDir = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir);
codeDir = fullfile(codeRootDir,'2pAnalysisCode');
analysisDir = fullfile(imagingRootDir,'analysis'); %directory of the analysis history summary file

addpath(genpath(codeRootDir)); % adds the paths of the code files to be executed. Useful on non-local computer such as the server

disp(['%%%%% CURRENT SESSION IS: ' animal ' ' session '. %%%%%']);

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

% ----- Logging of analysis ----- %
if ~exist(fullfile(imagingRootDir,'analysis'),'dir')
    mkdir(analysisDir);
end

try
    
    load(fullfile(analysisDir,'analysis.mat'),'analysis'); % load the analysis history summary file
    currentDate=datestr(now,'yyyy-mm-dd_HHMMss'); currAnalIdx=length(analysis)+1; % append data of the current analysis to the summary file
    
    disp('Loading 2P imaging data...');
    
    % organize suite2p data into npy files (ROI and spiking)
    if loadSavedData == false || nargout > 2
        disp('LOADING: npy file...'); tic;
        npy = twoP_importSuite2pNPY(animal,session);
        disp(['LOADED: The .npy file was loaded in ' num2str(toc) ' seconds.'])
        
        % IMPORTANT!!! smoothes inferred spikes with  a gaussian filter
        if doSmooth==true
            spks = smoothCol(npy.spks,2,filterWindow,'gauss');
        else
            spks = npy.spks;
        end
        
        % Creates stimulus-aligned data matrices and saves data.mat in the
        % same directory as npy files
        tic
        try
            data = twoP_alignDetectionTask(npy.ops, spks, npy.iscell, npy.redcell, npy.bin_MScan_filepath); % align suite2p data to sensory stimulus
            data.smoothed = doSmooth;
        catch ME
            data.ME = ME;
        end
        disp(['Trial alignment to stimulus onset was completed in ' num2str(toc) ' seconds.'])
        
        if exist('ME','var')
            disp(['Error detected: ' ME.identifier '. ' ME.message])
        else
            disp('2P DATA LOADED!');
        end
    else
        dataPath = fullfile(s2pDataDir,'data.mat');
        data = load(dataPath ,'data'); data = data.data;
    end
    
    idxRow = twoP_findRowAcquisitionRecord(animal,session);
    sFilenames = regexp(exps{idxRow,9},'\;','split','once');

    % Load behavior data
    disp('LOADING: behavior data...');
%     bhvDir = fullfile(bhvRootDir,animal,bhvSubDir);
    if length(sFilenames) > 1 % if the behavior session was broken into two due to bpod control computer crash
        [bhv,bhvFilePath] = twoP_RepairSession(animal,sFilenames,data.trialCodes);
    else
        [bhv,bhvFilePath] = twoP_loadBehaviorSession(animal,session);
    end
    disp(['LOADED: behavior data: ' bhvFilePath]);
    
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
    writetable(struct2table(analysis),fullfile(analysisDir, 'analysis.csv'));    
end