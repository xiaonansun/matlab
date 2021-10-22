function twoP_processAllSessions(varargin)
%%
if nargin > 0
    C=varargin{1};
end

docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';

% Specify session to load
exps = GetGoogleSpreadsheet(docid);
iStart = 2;

parfor i = iStart:size(exps,1)
    % for i = 95:114
    animal = exps{i,1};
    session = exps{i,6};
    
    try
        
        % ----- Specify session to load ----- %%
%         animal = 'CSP27';
%         session = '20200329';
        [data,SessionData]=twoP_loadImgBhvData(animal,session, true, 10, false);
        
        %
        D = twoP_combineStimAlignedData(data);
        D = twoP_adjustData(D,SessionData);
        
        % ----- Makes adjustments to the twoP data struct ----- %%
        
        cBhv = selectBehaviorTrials(SessionData, D.trialNumbers,animal, session); %% very important in matching trial indices
        
        Vc = rateDisc_getBhvRealignment(D.neural, cBhv, [], [], animal, session); % re-aligned imaging data to trial epoches
        
        twoP_compareStimTimes(D,cBhv);
    end
end
