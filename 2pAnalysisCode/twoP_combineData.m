% function twoP_combineData

% This function combines inferred spikes across sessions into a single
% cell by performing the following:
% 

exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1); colSession = exps(:,6);

S = twoP_settings;
nShuffle = S.nShuffle;
if S.isUnix == true; parpool('local',32); end
bhvRootDir = S.dir.bhvRootDir;
bhvSubDir = S.dir.bhvSubDir;
baseDir = S.dir.imagingRootDir; s2pDir = S.dir.imagingSubDir;
celltype = cell(size(exps,1),1);
trialTypes = {};
for i = 1:length(S.cellTypes)
    ctIdx = contains(exps(:,1),S.cellTypes{i});
    celltype(ctIdx) = S.cellTypes(i);
end
segIdx = S.segIdx;

D = cell(size(exps,1),10);
Vsub = cell(size(exps,1),8);
D{1,1} = 'animal'; D{1,2} = 'session'; D{1,3} = 'Last trialNumber';
D{1,4} = 'total trialNumbers'; D{1,5} = 'Number of trials in data.neural';
D{1,6} = 'Number of rewarded trials';
D{1,7} = 'Number of rewarded trials - corrected';

tic
parfor i = 2:size(exps,1) % This loop performs computations for individual sessions
% parfor i = 2:20
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
        rD{8} = Vc(data.idx_redcell,:,:);
        rD{9} = Vc(data.idx_notredcell,:,:);
        rD{10} = cBhv;
        D(i,:) = rD;
        disp(['Loaded ' animal ' ' session '.']);
        
        sBhv = twoP_bhvSubSelection(cBhv);
        
        disp('Computing mean PSTH...');
        rVsub = Vsub(i,:);
        rVsub{1} = zeros(length(data.idx_redcell),size(Vc,2),size(sBhv.sub.AllIdx,1)); % Pre-allocates memory
        rVsub{2} = zeros(length(data.idx_redcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        rVsub{3} = zeros(length(data.idx_notredcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        rVsub{4} = zeros(length(data.idx_notredcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        rVsub{5} = zeros(length(data.idx_redcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        rVsub{6} = zeros(length(data.idx_redcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        rVsub{7} = zeros(length(data.idx_notredcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        rVsub{8} = zeros(length(data.idx_notredcell),size(Vc,2),size(sBhv.sub.AllIdx,1));
        for iSub = 1:size(sBhv.sub.AllIdx,1)
            rVsub{1}(:,:,iSub) = squeeze(mean(Vc(data.idx_redcell,:,sBhv.sub.AllIdx(iSub,:)),3,'omitnan'));
            rVsub{2}(:,:,iSub) = squeeze(std(Vc(data.idx_redcell,:,sBhv.sub.AllIdx(iSub,:)),0,3,'omitnan'));
            rVsub{3}(:,:,iSub) = squeeze(mean(Vc(data.idx_notredcell,:,sBhv.sub.AllIdx(iSub,:)),3,'omitnan'));
            rVsub{4}(:,:,iSub) = squeeze(std(Vc(data.idx_notredcell,:,sBhv.sub.AllIdx(iSub,:)),0,3,'omitnan'));
            tempR = zeros(length(data.idx_redcell),size(Vc,2),nShuffle);
            tempNR = zeros(length(data.idx_notredcell),size(Vc,2),nShuffle);
            for iShuf = 1:nShuffle
                shufIdx = randperm(size(sBhv.sub.AllIdx,2));
                tempR(:,:,iShuf) = squeeze(mean(Vc(data.idx_redcell,:,sBhv.sub.AllIdx(iSub,shufIdx)),3,'omitnan'));
                tempNR(:,:,iShuf) = squeeze(mean(Vc(data.idx_notredcell,:,sBhv.sub.AllIdx(iSub,shufIdx)),3,'omitnan'));
            end
            rVsub{5}(:,:,iSub) = mean(tempR,3,'omitnan');
            rVsub{6}(:,:,iSub) = std(tempR,0,3,'omitnan');
            rVsub{7}(:,:,iSub) = mean(tempNR,3,'omitnan');
            rVsub{8}(:,:,iSub) = std(tempNR,0,3,'omitnan');
        end
        Vsub(i,:) = rVsub;
        trialTypes(i,:) = sBhv.sub.names;
        disp('Computing mean PSTH... DONE!')
    end
end
Vsub = [Vsub celltype];
% trialTypes = sBhv.sub.names;

disp(['All sessions combined in ' num2str(toc) ' seconds.']);

% T = cell2table(D);
% writetable(T,fullfile(S.dir.imagingRootDir,'analysis','trialNumberComparison.csv'));

% save(fullfile(S.dir.imagingRootDir,'analysis','aligned_combined.mat'),'D','-nocompression','-v7.3');
tic
save(fullfile(S.dir.imagingRootDir,'analysis','all_psth.mat'),'Vsub','trialTypes','-nocompression','-v7.3');
disp(['Data saved in' num2str(toc) ' seconds.']);