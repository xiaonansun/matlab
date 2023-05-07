function twoP_processAllSessions(varargin)
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
        animal = 'Plex51'; session = '200320'; % <---- modify this line to load one session
        procSingleSession = 1;
    else
        animal = colAnimal{i};
        session = colSession{i};
    end
    
    try
        % NOTE: When loading isredcell, this must be loaded from the npy
        % files, not from any previously-saved .mat files. If anything
        % changed with cell curation, loading .mat would not be able to
        % load the updated neuron definitions.
        
        [data,SessionData]=twoP_loadImgBhvData(animal,session, true, 10, false);
%         data = twoP_zscoreTrializeContinuousSDU(animal, session, true);

        D = twoP_combineStimAlignedData(data); % Combines fields of data struct that are organized in cells (this occurs when imaging was interrupted by an MScan crash)
        D = twoP_adjustData(D,SessionData);
        
        cBhv = selectBehaviorTrials(SessionData, D.trialNumbers,animal, session); %% very important in matching trial indices
        
        rateDisc_getBhvRealignment(D.neural, cBhv, [], [], animal, session); % Re-aligns imaging data to epoches
        rateDisc_getBhvRealignment(D.sdu, cBhv, [], [], animal, session, 'sduVc'); % Re-aligns z-scored imaging data to epoches
    
        % REQUEST: A histogram should be added to the right side of the
        % scatter plot. This is necessary since there may be a few black
        % points buried within the red and blue points that may be subtle.
        % Having black data points distributed within the red and blue
        % data indicates some twoP data not matched to the bpod data.
        twoP_compareStimTimes(D,cBhv);
    catch ME
        disp(ME.message);
    end
    
    if ~exist('procSingleSession','var') || procSingleSession == 0
        close all;
    end
end