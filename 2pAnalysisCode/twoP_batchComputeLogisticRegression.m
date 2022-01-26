function twoP_batchComputeLogisticRegression(varargin)
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

for i = iStart:size(exps,1)
    %% execute this cell to load one session. Just modify the animal and session variables as needed
    
    try

        S = twoP_settings;
        imagingRootDir = S.dir.imagingRootDir;
        imagingSubDir = S.dir.imagingSubDir;

        animal = colAnimal{i};
        session = colSession{i};

        twoP_computeLogisticRegression(animal,session,50);
        
%         twoP_aucAnalysisNew(Vc, cBhv, 1, 1, animal, session);
        
    end
%     if ~exist('procSingleSession','var') || procSingleSession == 0
%         close all;
%     end
end
