function rateDisc_twoP_locaChoiceDecode(animal)

%% some basic variables
stepSize = 15;
regType = 'lasso'; %lasso or ridge
learnType = 'leastsquares';
sRate = 30;

if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for global dimensions

%% load data and determine training ranges
load([bPath 'mask.mat'],'mask');
load([bPath 'locAC.mat'],'A', 'C', 'nanIdx');
load([bPath 'trialInfo.mat'],'bhvTrials','recs','trialCnt');
load([bPath 'opts.mat'],'opts');
load([bPath animal '_blockBhv.mat'],'bhv');

frameCnt = length(nanIdx) / sum(trialCnt); %number of frames per trial
A = arrayShrink(A, mask, 'merge');

[dataOverview, ~, ~, ~, segIdx, segLabels, ~, ~, trainDates] = rateDiscRecordings;
recIdx = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs); %for each recording, determine to which part of the training it belongs
if any(recIdx == 0)
    fprintf('%d/%d recordings are unassigned\n',sum(recIdx==0),length(recIdx))
end
segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment

load('allenDorsalMapSM.mat')
for iAreas = 1 : length(dorsalMaps.edgeOutlineSplit)
    dorsalMaps.edgeOutlineSplit{iAreas}(:,1) = dorsalMaps.edgeOutlineSplit{iAreas}(:,1)-yRange(1);
    dorsalMaps.edgeOutlineSplit{iAreas}(:,2) = dorsalMaps.edgeOutlineSplit{iAreas}(:,2)-xRange(1);
end

%% isolate data subset for each training set, align to trial events and run decoder
allTrials = [0 cumsum(trialCnt)];
for iPhase = unique(recIdx(recIdx ~= 0)) % don't use 0s
    
    tIdx = allTrials(find(recIdx == iPhase,1,'first') : find(recIdx == iPhase,1,'last')+1); %trial idx
    cBhv = selectBehaviorTrials(bhv, tIdx(1)+1 : tIdx(end)); %select trials from current learning phase
    
    Vc = NaN(size(C,1), frameCnt*length(cBhv.Rewarded), 'single');
    firstFrame = tIdx(1)*frameCnt+1;
    for iRecs = 1 : length(tIdx)-1 % loop through recordings to collect data and normalize
        
        fIdx = tIdx(iRecs)*frameCnt+1 : tIdx(iRecs+1)*frameCnt; %frame idx
        cData = C(:, sum(~nanIdx(1:fIdx(1)-1)) + (1:sum(~nanIdx(fIdx)))); %data from current session
%         cData = zscore(cData, [], 2); %z-score all components
        targIdx = fIdx - firstFrame + 1; %index for where to place data in Vc array
        Vc(:,targIdx(~nanIdx(fIdx))) = cData; %fill non-frames from C into current selection
        
    end
    
    Vc = reshape(Vc,size(Vc,1), frameCnt, []);
    Vc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts); %aligned to different trial episodes
    % Vc: neural data neurons x frame x trial
    
[cvAcc, bMaps, trialCnt] = rateDisc_logDecoder(Vc, [], cBhv, 400, 0, regType); 
% [cvAcc, bMaps, trialCnt] = rateDisc_logDecoder(data, U, bhv, useTrials, targMod, regType, stepSize, decType)
% useTrials: 400 trials is default, can be adjusted, fewer trials will reduce the accuracy of the decoder
% stepSize: will downsample data
% decType: the type of decoder, default is allchoice

end