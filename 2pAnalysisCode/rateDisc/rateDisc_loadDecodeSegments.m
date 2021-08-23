% function rateDisc_loadDecodeSegments

%% get some general information
% animals = {'mSM63' 'mSM64' 'mSM65' 'mSM66'};
animals = {'mSM66'};
if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end

[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
load('allenDorsalMapSM.mat'); %get area outlines and match mask size
for iAreas = 1 : length(dorsalMaps.edgeOutlineSplit)
    dorsalMaps.edgeOutlineSplit{iAreas}(:,1) = dorsalMaps.edgeOutlineSplit{iAreas}(:,1);
    dorsalMaps.edgeOutlineSplit{iAreas}(:,2) = dorsalMaps.edgeOutlineSplit{iAreas}(:,2);
end
decType = 'stim'; %determines which trials were used for decoder. use 'all','correct','error','choice' or 'stim'
cMod = 'all'; %target modality. use 'all','audio','tactile' or 'audiotactile'.

betaMaps = cell(1, length(animals));
decAcc = cell(1, length(animals));
trialCnts = cell(1, length(animals));

%% loop through animals
for iAnimal = 1 : length(animals)
    animal = animals{iAnimal};
    
%     bPath = [cPath animal filesep 'blockData' filesep 'nmfDecSegments' filesep]; % path for global dimensions
    bPath = [cPath animal filesep 'blockData' filesep 'decodeSegments' filesep]; % path for global dimensions
    load([cPath animal filesep 'blockData' filesep 'mask.mat'])
    load([cPath animal filesep 'blockData' filesep 'trialInfo.mat'])
    allTrials = [0 cumsum(trialCnt)]; %max(recIdex) per sessions
    recIdx = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs); %for each recording, determine to which part of the training it belongs
    
    ranges = []; %define range for detection training and discrimination. this depends on the chosen modality
    if strcmpi(cMod, 'audio')
        ranges = {1:2 3:4 5}; %training range for audio data. Entries define range for early/late detection and discrmination
    elseif strcmpi(cMod, 'tactile')
        ranges = {6:7 8:9 [6 10]}; %training range for tactile data. Entries define range for early/late detection and novice/expert discrmination
    end
    
    %% load segments
    files = dir([bPath cMod '_' decType '_seg*']);
    cFile = matfile([bPath files(1).name]);
    betaMaps{iAnimal} = NaN([size(cFile.bMaps), max(recIdx)], 'single'); %running average for betamaps at different training stages
    decAcc{iAnimal} = NaN([size(cFile.cvAcc,2) length(files)], 'single'); %decoder accuracy
    trialCnts{iAnimal} = NaN([size(cFile.cvAcc,2) length(files)], 'single'); %trialcounts for different bins
    if ~isempty(ranges)
        rangeMaps{iAnimal} = NaN([size(cFile.bMaps), length(ranges)], 'single'); %running average for betamaps at different training stages
    else
        rangeMaps = NaN;
    end
    clear cFile
    
    recCnt = zeros(1,max(recIdx)); %counter for running average in different recIdx entries
    rangeCnt = zeros(1,length(ranges)); %counter for running average in different training episodes
    recIDs{iAnimal} = NaN(1,length(files)); %identifier for which segments were added to which recIdx entry
    rangeIDs{iAnimal} = NaN(1,length(files)); %identifier for which segments were added to training ranges
    allIdx = NaN(2, length(files), 'single');
    
    for iSegs = 1 : length(files)
        % for iSegs = 1 : 5
        load([bPath files(iSegs).name], 'bMaps', 'cvAcc','cIdx','trialCntWin');
        allIdx(:,iSegs) = cIdx;
        
        % check where this segment belongs
        rIdx = max(recIdx(allTrials < cIdx(1) + round(diff(cIdx)/2)));
        betaMaps{iAnimal}(:,:,rIdx) = nansum(cat(3,betaMaps{iAnimal}(:,:,rIdx)*recCnt(rIdx), bMaps),3) ./ (recCnt(rIdx) + 1); %running average
        recCnt(rIdx) = recCnt(rIdx) + 1;
        recIDs{iAnimal}(iSegs) = rIdx;
        
        % check for range
        for iRange = 1:length(ranges)
            if rIdx >= ranges{iRange}(1) && rIdx <= ranges{iRange}(end) %part of current range
                rangeMaps{iAnimal}(:,:,iRange) = nansum(cat(3,rangeMaps{iAnimal}(:,:,iRange)*rangeCnt(iRange), bMaps),3) ./ (rangeCnt(iRange)+1); %running average
                rangeCnt(iRange) = rangeCnt(iRange) + 1;
                rangeIDs{iAnimal}(iSegs) = iRange;
            end
        end
        
        % keep classifier accuracy and trialcounts
        decAcc{iAnimal}(:,iSegs) = cvAcc;
        trialCnts{iAnimal}(:,iSegs) = trialCntWin;
    end
    a = arrayShrink(betaMaps{iAnimal},mask,'split');
    betaMaps{iAnimal} = NaN([size(allenMask),size(a,3),size(a,4)],'single');
    betaMaps{iAnimal}(yRange,xRange,:,:) = a; clear a
    betaMaps{iAnimal} = arrayShrink(betaMaps{iAnimal},allenMask);
    
    if iscell(rangeMaps)
        a = arrayShrink(rangeMaps{iAnimal},mask,'split');
        rangeMaps{iAnimal} = NaN([size(allenMask),size(a,3),size(a,4)],'single');
        rangeMaps{iAnimal}(yRange,xRange,:,:) = a; clear a
        rangeMaps{iAnimal} = arrayShrink(rangeMaps{iAnimal},allenMask);
    end
end

%%
% cMaps = mean(betaMaps(:,:,~isnan(betaMaps(1,90,:))),3);
% cMaps = betaMaps(:,:,2);
targRec = 1:10; %auditory discrimination
cMaps = [];
cInd = false(size(trialCnts{1},1),1);
for iAnimal = 1 : length(animals)
    recCnt = zeros(1,length(targRec));
    cData = betaMaps{iAnimal}(:, :, targRec);
    for iRec = 1 : length(targRec)
        recCnt(iRec) = sum(ismember(recIDs{iAnimal},targRec(iRec))); %check how many segments are part of current range
        cData(:,:,iRec) = cData(:,:,iRec) .* recCnt(iRec);
    end
    cData = nansum(cData,3) ./ sum(recCnt); %count-corrected average
    cInd = cInd | mean(trialCnts{iAnimal}(:, ismember(recIDs{iAnimal},targRec)),2) < 200; %dont use bins with less than 250 trials
    cData(:,cInd) = NaN;
    cMaps = cat(3,cMaps,cData);
end
clear cData
cMaps = nanmean(cMaps,3);
cMaps(cMaps == 0) = NaN;

% cMaps = cat(3,betaMaps{:});
% b = mean(cat(2,trialCnts{:}),2) < 250;
% cMaps = nanmean(cMaps,3);
% cMaps(:,b) = NaN;
% clear nanIdx
% nanIdx(1,:) = [0 find(diff(~(isnan(cMaps(1,:)))) == 1)]+1;
% if isnan(cMaps(1,end))
%     nanIdx(2,:) = find(diff(~(isnan(cMaps(1,:)))) == -1);
% else
%     nanIdx(2,:) = [find(diff(~(isnan(cMaps(1,:)))) == -1) size(cMaps,2)];
% end
% for x = 1 : size(nanIdx,2)
%     cMaps(:, nanIdx(1,x):nanIdx(2,x)) = smoothCol(cMaps(:, nanIdx(1,x):nanIdx(2,x)),3,'gauss',[],2);
% end

cMaps = arrayShrink(cMaps,allenMask,'split');
temp = cMaps(:,1:(size(cMaps,2)/2),:) - cMaps(:,end : -1 : (size(cMaps,2)/2) + 1,:);


% end