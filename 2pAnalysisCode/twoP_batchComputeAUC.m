function twoP_batchComputeAUC(varargin)
% This is a function 

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
    
%     if ~exist('N','var') || N == 0
%         animal = 'Plex50'; session = '200322'; % <---- modify this line to load one session
%         procSingleSession = 1;
%     else
%         animal = colAnimal{i};
%         session = colSession{i};
%     end
    
    try
        
        %         if ~exist('N','var') || N == 0
        %         animal = 'CSP27'; session = '20200319';
        %         end

        S = twoP_settings;
        imagingRootDir = S.dir.imagingRootDir;
        imagingSubDir = S.dir.imagingSubDir;

        animal = colAnimal{i};
        session = colSession{i};

        % Load event-aligned imaging data, behavior data, cell-type ID data
        Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mat'),'Vc'); Vc = Vc.Vc;
        cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
%         idxCell = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'iscell.npy'));
%         idxRed = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'redcell.npy'));
%         idxRed = logical(idxRed(logical(idxCell(:,1))));
        
        % Define the epoches
%         idxEpochAll = [idxEpochInit; idxEpochStimEarly; idxEpochStimLate; idxEpochDelay; idxEpochResponseEarly; idxEpochResponseLate];

        % compute trials binned by epochs
        eVc = twoP_epochTrialMean(Vc,S.allEpoches.idx);
        
        % compute and save epoch-binned and unbinned AUC
        twoP_aucAnalysisNew(eVc, cBhv, 1, 1, animal, session);
        
        twoP_aucAnalysisNew(Vc, cBhv, 1, 1, animal, session);
        
    end
%     if ~exist('procSingleSession','var') || procSingleSession == 0
%         close all;
%     end
end
