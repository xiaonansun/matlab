function rateDisc_nmfChoiceDecode(animal, stepWin, decType, cMod, trialsPerBin)

%% some basic variables
tic;
if ~exist('decType','var') || isempty(decType)
    decType = 'all'; %all trials is default
end
if ~exist('cMod','var') || isempty(cMod)
    cMod = 2; %audio is default
end

if cMod == 0
    modTitle = 'all';
elseif cMod == 2
    modTitle = 'audio';
elseif cMod == 4
    modTitle = 'tactile';
elseif cMod == 6
    modTitle = 'audiotactile';
else
    error('unknown modality');
end
    
regType = 'lasso'; %lasso or ridge
sRate = 30;

if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for global dimensions

%% load data and determine training ranges
load([bPath 'mask.mat'],'mask','xRange','yRange');
load([bPath 'trialInfo.mat'],'bhvTrials','recs','trialCnt');
load([bPath 'opts.mat'],'opts');
load([bPath animal '_blockBhv.mat'],'bhv');

[dataOverview, ~, ~, ~, segIdx, ~, ~, ~, trainDates] = rateDiscRecordings;
recIdx = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs); %for each recording, determine to which part of the training it belongs
if any(recIdx == 0)
    fprintf('%d/%d recordings are unassigned\n',sum(recIdx==0),length(recIdx))
end
segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment

load('allenDorsalMapSM.mat'); %get area outlines and match mask size
for iAreas = 1 : length(dorsalMaps.edgeOutlineSplit)
    dorsalMaps.edgeOutlineSplit{iAreas}(:,1) = dorsalMaps.edgeOutlineSplit{iAreas}(:,1)-yRange(1);
    dorsalMaps.edgeOutlineSplit{iAreas}(:,2) = dorsalMaps.edgeOutlineSplit{iAreas}(:,2)-xRange(1);
end

load([bPath 'locAC.mat'],'A','nanIdx');

%% isolate data subset for each training set, align to trial events and run decoder
allTrials = [0 cumsum(trialCnt)]; %trials per sessions
frameCnt = size(nanIdx,2) / allTrials(end); %frames per trial

% for each window, find the required amount of cases when combining choice, stimulus and correct/error
leftIdx = (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left (choice)
corrIdx = bhv.Rewarded; %rewarded trials
if cMod == 0
    modIdx = true(1,length(bhv.StimType)); %use all trials
else
    modIdx =  bhv.StimType == cMod; %select modality
end

% get combined cases to determine window size
binFac = 1.5; %by default, use three times more trials as required to get trials in the delay
clear tCombs
if strcmpi(decType, 'all')
    tCombs{1} = find(leftIdx & modIdx);
    tCombs{2} = find(~leftIdx & modIdx);
elseif strcmpi(decType, 'correct')
    tCombs{1} = find(leftIdx & modIdx & corrIdx);
    tCombs{2} = find(~leftIdx & modIdx & corrIdx);
elseif strcmpi(decType, 'error')
    tCombs{1} = find(leftIdx & modIdx & ~corrIdx);
    tCombs{2} = find(~leftIdx & modIdx & ~corrIdx);
elseif strcmpi(decType, 'choice') || strcmpi(decType, 'stim')
    tCombs{1} = find(corrIdx & leftIdx & modIdx);
    tCombs{2} = find(corrIdx & ~leftIdx & modIdx);
    tCombs{3} = find(~corrIdx & ~leftIdx & modIdx);
    tCombs{4} = find(~corrIdx & leftIdx & modIdx);
    binFac = binFac/2; %this needs to be adjusted to account for using 4 cases instead of 2
end

% get index for all steps until current step
lastTrial = 0; breaker = false;
for iSteps = 1 : stepWin
    cIdx = [lastTrial lastTrial];
    for x = 1 : length(tCombs)
       temp = tCombs{x}(tCombs{x} > lastTrial); %find the case that needs the most trials to reach the required case count
       if length(temp) > ceil(trialsPerBin*binFac)
           y = temp(ceil(trialsPerBin*binFac));
       elseif ~isempty(temp)
           y = temp(end);
           breaker = true;
       end
       if y > cIdx(2)
           cIdx(2) = y; %use selected case to find the end of the trialindex
       end
    end
    lastTrial = lastTrial + round(diff(cIdx) / 5); %move with 80% overlap
    if breaker; break; end
end
fprintf([animal ' - Current step: %d\n'], stepWin)
fprintf('Behavior loaded\n');
toc;

if stepWin <= iSteps
    cBhv = selectBehaviorTrials(bhv, cIdx(1)+1:cIdx(2)); %behavior for current window
    
    
    fIdx = [cIdx(1)*frameCnt+1 (cIdx(2))*frameCnt]; %frame idx
    sourceIdx = (fIdx(1) - sum(nanIdx(1:fIdx(1)-1))) + (0 : sum(~nanIdx(fIdx(1):fIdx(2)))-1); %index in C

    cFile = matfile([bPath 'locAC.mat']); %check wV file to load subset of components
    Vc = NaN(size(A,2), frameCnt*diff(cIdx),'single'); %pre-allocate dataset    
    Vc(:, ~nanIdx(fIdx(1) : fIdx(2))) = cFile.C(1:size(A,2), sourceIdx); %data for current sessions
    Vc = reshape(Vc,size(Vc,1),frameCnt,[]);
    Vc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts); %aligned to different trial episodes
    [cvAcc, bMaps, trialCntWin] = rateDisc_logDecoder(Vc, A, cBhv, trialsPerBin, cMod, regType, decType);
    fprintf('Decoder finished\n')
    toc;

    if ~exist([bPath 'nmfDecSegments'], 'dir'); mkdir([bPath 'nmfDecSegments']); end
    save([bPath 'nmfDecSegments' filesep modTitle '_' decType '_seg_' num2str(stepWin, '%04i')], 'cvAcc', 'bMaps', 'trialCntWin', 'cIdx', 'allTrials', 'recs','-v7.3')
    fprintf('Results saved\n')
    toc;
end

% allTrials = [0 cumsum(trialCnt)];
% for iPhase = unique(recIdx(recIdx ~= 0)) % don't use 0s
%     
%     tIdx = allTrials(find(recIdx == iPhase,1,'first') : find(recIdx == iPhase,1,'last')+1); %trial idx
%     cBhv = selectBehaviorTrials(bhv, tIdx(1)+1 : tIdx(end)); %select trials from current learning phase
%     
%     Vc = NaN(size(wV,1), frameCnt*length(cBhv.Rewarded), 'single');
%     Cnt = 0;
%     for iRecs = 1 : length(tIdx)-1 % loop through recordings to collect data and normalize
%         fIdx = tIdx(iRecs)*frameCnt+1 : tIdx(iRecs+1)*frameCnt; %frame idx
%         cData = wV(:, fIdx); %data from current session
% %         cData(:,~isnan(cData(1,:))) = zscore(cData(:,~isnan(cData(1,:))), [], 2); %z-score all components
% %         cData = bsxfun(@rdivide,cData,nanstd(cData, 0, 2)); %z-score all components
%         Vc(:, Cnt + (1:length(fIdx))) = cData;
%         Cnt = Cnt + length(fIdx);
%     end
%     
%     Vc = reshape(Vc,size(Vc,1), frameCnt, []);
%     Vc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts); %aligned to different trial episodes
%     [cvAcc, bMaps, trialCnt] = rateDisc_logDecoder(Vc, U, cBhv, 500, 0, regType);
% end

