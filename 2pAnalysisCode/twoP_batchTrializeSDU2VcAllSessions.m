function twoP_batchTrializeSDU2VcAllSessions

%% One way to trialize the new data: take the frame number of stimulus onset (data.stimSamplesOrig) to extract single trials from the continuous sdu trace

S = twoP_settings;
exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1);
colSession = exps(:,6);
imagingRootDir = S.dir.imagingRootDir;
% imagingSubDir = S.dir.imagingSubDir;
error_log = cell(size(exps,1),1);

%%
parfor i = 1:size(exps,1)
    %% i=141
    animal = colAnimal{i};
    session = colSession{i};
%     data_path = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'data.mat');
%     sdu_Vc_savepath = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'sduVc.mat');
%     bhv_path = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat');
    try
%     data = load(data_path); data = data.data;
%     cBhv = load(bhv_path); cBhv = cBhv.cBhv;
        
    twoP_zscoreTrializeContinuousSDU(animal, session, false);

    [data,SessionData]=twoP_loadImgBhvData(animal,session, true, 10, true);
    
    D = twoP_combineStimAlignedData(data); % Combines fields of data struct that are organized in cells (this occurs when imaging was interrupted by an MScan crash) into matrices
    D = twoP_adjustData(D,SessionData); 
    
    cBhv = selectBehaviorTrials(SessionData, D.trialNumbers,animal, session); %% very important in matching trial indices
        
    rateDisc_getBhvRealignment(D.sdu, cBhv, [], [], animal, session, 'sduVc');
%     save(sdu_Vc_savepath,'sduVc'); disp(['Saved sduVc as: ' sdu_Vc_savepath]);
    catch ME
        disp(ME.message)
        error_log{i}= ME.message;
    end
    
end

st = dbstack; % get the file name of the current function
if ~isempty(st)
namestr = st.name; clear st;
save(fullfile(imagingRootDir,[namestr '_error_log.mat']),'error_log');
end