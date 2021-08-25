% function twoP_combineData

exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1); colSession = exps(:,6);

S = twoP_settings;
bhvRootDir = S.dir.bhvRootDir;
bhvSubDir = S.dir.bhvSubDir;
baseDir = S.dir.imagingRootDir; s2pDir = S.dir.imagingSubDir;

segIdx = S.segIdx;

D = cell(size(exps,1),10);
D{1,1} = 'animal'; D{1,2} = 'session'; D{1,3} = 'Last trialNumber';
D{1,4} = 'total trialNumbers'; D{1,5} = 'Number of trials in data.neural';
D{1,6} = 'Number of rewarded trials';
D{1,7} = 'Number of rewarded trials - corrected';

tic
parfor i = 2:size(exps,1)
% parfor i = 2:10
    try
        animal = colAnimal{i};
        session = colSession{i};
        %         [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10, false);
        dataPath = fullfile(baseDir,animal,'imaging',session,s2pDir,'data.mat');
        bhvDir = fullfile(bhvRootDir,animal,bhvSubDir);
        disp(['Loading ' animal ' ' session '...']);
        
        data = load(dataPath ,'data'); data = data.data;
        
        bhv = twoP_loadBehaviorSession(animal,session,bhvDir);
        rD = D(i,:);
        rD{1} = animal; rD{2} = session; rD{3} =  max(data.trialNumbers);
        rD{4} = length(data.trialNumbers); rD{5} = size(data.neural,3);
        rD{6} = length(bhv.Rewarded);
        
        data = twoP_adjustData(data,bhv);
        
        opts = struct;
        opts.preStim = data.trialStimFrame*data.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
        opts.frameRate = 1000/data.msPerFrame; % Frame rate of imaging
        sRate = opts.frameRate;
        segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment
        
        cBhv = selectBehaviorTrials(bhv, data.trialNumbers); %% Match trial indices - IMPORTANT!
        rD{7} = length(cBhv.Rewarded);
        
        
        Vc = rateDisc_getBhvRealignment(data.neural, cBhv, segFrames, opts); % re-aligned imaging data to trial epoches
        %         D{i,1} = Vc(data.idx_redcell,:,:);
        rD{8} = Vc(data.idx_redcell,:,:);
        rD{9} = Vc(data.idx_notredcell,:,:);
        rD{10} = cBhv;
        
        D(i,:) = rD;
        
        disp(['Loaded ' animal ' ' session '.']);
        
    end
end

disp(['All sessions loaded in ' num2str(toc) ' seconds.']);

% T = cell2table(D);
% writetable(T,fullfile(S.dir.imagingRootDir,'analysis','trialNumberComparison.csv'));
save(fullfile(S.dir.imagingRootDir,'analysis','aligned_combined.mat'),'D','-nocompression','-v7.3');