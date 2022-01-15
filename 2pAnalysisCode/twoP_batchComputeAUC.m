function twoP_batchComputeAUC(varargin)
% This is a batch function that processes suite2p output (i.e. .npy
% files) as well as the corresponding Bpod behavior data into trialized
% format.
% Included functions are:
% (1) twoP_loadImgBhvData.m: Converts npy files into trialized imaging data
% where the dimensions of the output (data.neural) is neurons x frame x
% trial.
% (2) twoP_adjustData.m: Makes minor adjustments to the data struct:
% Removes extra data.trialNumbers if length(data.trialNumbers) >
% size(data.neural,3). Removes extra data.trialNumbers if length(data.trialNumbers) >
% length(SessionData.Rewarded). Removes extra data.neural trials if
% size(data.neural,3) > length(data.trialNumbers).
% (3) twoP_combineStimAlignedData.m: Some sessions are interrupted
% by MScan computer crashes, essentially splitting the MScan .bin file into
% two. The previous function, twoP_loadImgBhvData.m, adjusts for this
% possibility by splitting the overall session data into multiple cells
% that correspond to the number of MScan .bin files that may exist (usually
% two files). This function then checks the trial numbers, deletes
% incomplete trials, and combines multiple cells of data into a single
% matrix.
% (4) selectBehaviorTrials.m: an additional step in matching trials between
% the behavior data (Bpod output) and twoP data (data struct).
% (5) rateDisc_getBhvRealignment.m: Further aligns the data.neural struct
% to epoches of the trial: handle initation, simulus onset, delay epoch
% onset, and response onset.
% (6) twoP_compareStimTimes.m: plots the inter-trial interval (the
% duration between the stimulus onset of consecutive trials) of the MScan
% timestamp versus the Bpod timestamp. If all values fall along x=y, then
% then trial numbers are matched between behavior and twoP data.

%%

if nargin > 0
    C=varargin{1};
end
N = nargin;

docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';

% Specify session to load
exps = GetGoogleSpreadsheet(docid);
colAnimal = exps(:,1);
colSession = exps(:,6);

iStart = 2;

parfor i = iStart:size(exps,1)
    %% execute this cell to load one session. Just modify the animal and session variables as needed
    
    if ~exist('N','var') || N == 0
        animal = 'Plex50'; session = '200322'; % <---- modify this line to load one session
        procSingleSession = 1;
    else
        animal = colAnimal{i};
        session = colSession{i};
    end
    
    try
        
        %         if ~exist('N','var') || N == 0
        %         animal = 'CSP27'; session = '20200319';
        %         end
        
        S = twoP_settings;
        imagingRootDir = S.dir.imagingRootDir;
        imagingSubDir = S.dir.imagingSubDir;
        
        sPerFrame = S.msPerFrame/1000;
        
        % Load event-aligned imaging data, behavior data, cell-type ID data
        Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mat'),'Vc'); Vc = Vc.Vc;
        cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
        idxCell = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'iscell.npy'));
        idxRed = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'redcell.npy'));
        idxRed = logical(idxRed(logical(idxCell(:,1))));
        
        % Define the epoches
        idxEpochInit = [S.segFrames(1) S.segFrames(2)];
        idxEpochStim = [S.segFrames(2)+1 S.segFrames(3)];
        idxEpochStimEarly = [S.segFrames(2)+1 floor(mean(idxEpochStim))];
        idxEpochStimLate = [ceil(mean(idxEpochStim)) S.segFrames(3)];
        idxEpochDelay = [S.segFrames(3)+1 S.segFrames(4)];
        idxEpochResponse = [S.segFrames(4)+1 S.segFrames(5)];
        idxEpochResponseEarly = [S.segFrames(4)+1 floor(mean(idxEpochResponse))];
        idxEpochResponseLate = [ceil(mean(idxEpochResponse)) S.segFrames(5)];
        idxEpochAll = [idxEpochInit; idxEpochStimEarly; idxEpochStimLate; idxEpochDelay; idxEpochResponseEarly; idxEpochResponseLate];

        % compute trials binned by epochs
        eVc = twoP_epochTrialMean(Vc,idxEpochAll);
        
        % compute and save epoch-binned and unbinned AUC
        [eAUC.all, eAUC.shuffle, eAUC.sAUC] = twoP_aucAnalysisNew(eVc, cBhv, 1, 1);
        parsave(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'eAUC.mat'),'eAUC');
        
        [AUC.all, AUC.shuffle, AUC.sAUC] = twoP_aucAnalysisNew(Vc, cBhv, 1, 1);
        parsave(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'AUC.mat'),'eAUC');
        
    end
    if ~exist('procSingleSession','var') || procSingleSession == 0
        close all;
    end
end
