% function twoP_combineData

exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1); colSession = exps(:,6);

S = twoP_settings;
baseDir = S.dir.imagingRootDir; s2pDir = S.dir.imagingSubDir;
bhvDir = fullfile(S.dir.bhvRootDir,animal,S.dir.bhvSubDir);
segIdx = S.segIdx;

D = cell(size(exps,1),2);

tic
parfor i = 2:size(exps,1)
    try
        animal = colAnimal{i};
        session = colSession{i};
        %         [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10, false);
        dataPath = fullfile(baseDir,animal,'imaging',session,s2pDir,'data.mat');
        
        disp(['Loading ' animal ' ' session '...']);
        
        data = load(dataPath ,'data'); data = data.data;
        [bhv,behaviorFilePath] = twoP_loadBehaviorSession(animal,session,bhvDir)
        opts = struct;
        opts.preStim = data.trialStimFrame*data.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
        opts.frameRate = 1000/data.msPerFrame; % Frame rate of imaging
        sRate = opts.frameRate;
        segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment
        
        cBhv = selectBehaviorTrials(SessionData, data.trialNumbers); %% Match trial indices - IMPORTANT!
        Vc = rateDisc_getBhvRealignment(data.neural, cBhv, segFrames, opts); % re-aligned imaging data to trial epoches
        D{i,1} = Vc(data.idx_redcell,:,:);

        %         D{i} = data.neural(1,1,1);
        disp(['Loaded ' animal ' ' session '.']);
        
    end
end

disp(['All sessions loaded in ' num2str(toc) ' seconds.']);