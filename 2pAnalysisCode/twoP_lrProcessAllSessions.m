function twoP_lrProcessAllSessions(varargin)

if nargin > 0
    C=varargin{1};
end

docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';

% Specify session to load
exps = GetGoogleSpreadsheet(docid);
iStart = 2;

for i = iStart:size(exps,1)
% for i = 95:114
    animal = exps{i,1};
    session = exps{i,6};
    
    try
%         baseFileName = [animal '_' session];
        [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10, false);
        
%         [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10);
        
        % ----- makes adjustments to data struct ----- %%
        data = twoP_adjustData(data,SessionData);
        
        % ----- Define event-aligned matrices ----- %%
        [lick,data.lickWinIdx,data.lickWinMs,data.dataLick, data.dataLickTrialNumbers]=twoP_alignToLick(data, SessionData);
        
        lr = twoP_logisticRegression(animal, session, data, SessionData, true, true);
%         lr = twoP_logisticRegression(animal, session, data, SessionData, true);
    end
    
    clear data lr; close all;
end
