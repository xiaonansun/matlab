% function [recV, semV, allBeta, recLabels, dataPath, allmodIdx, allSideIdx, alignIdx, allRecIdx, baseLength, frames, stimTimes, cellCnt] = twoP_motorReconstruct(cMod)
% code to create reconstructed 2pimaging data based on motor regressors. cMod
% informs over which animals to use ('Visual', 'Audio' or 'All'). Last
% dimension of recV is different modalities, 1 is corr. vision, 2 is corr. audio, 
% 3 is all corr. trials. 4:6 is the same thing but for all trials.

%% select data sets
dataOverview = twoP_delayDecRecordings; %get overview for all recordings
[~, motorLabels] = delayDecRecordings; %get index for motor regressors
opMotorLabels = {'lLick' 'rLick' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'}; %operant motor regressors

animals = dataOverview(:,1);
Cnt = 0;
baseLength = inf;
postLength = inf;
rng(1) % For reproducibility

for iAnimals = 1:length(animals)
    
    %% load data
    fPath = dataOverview{iAnimals,5};
    Cnt = Cnt +1;
    fPath = dataOverview{iAnimals,5};
    dataPath{Cnt} = fPath; %store current data path
    load([fPath 'data.mat']);
    
    % load behavior and get modality indices
    bhvFile = strsplit(fPath,filesep);
    bhvFile = dir([fPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
    load([fPath bhvFile.name]);
    
    bTrials = data.trialNumbers;
    trials = bTrials;
    bTrials(~ismember(data.trialNumbers,data.bhvTrials)) = []; %don't use trials that have problems with trial onset times
    bTrials(SessionData.DidNotChoose(bTrials) | SessionData.DidNotLever(bTrials) | ~SessionData.Assisted(bTrials)) = []; %don't use unperformed/assisted trials
    
    data.dFOF(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    data.DS(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    data.analog(:,:,~ismember(data.trialNumbers,bTrials)) = [];
  
    clear modIdx
    sucInd = SessionData.Rewarded(bTrials) & SessionData.Assisted(bTrials); %find succesful unisensory trials
    modIdx{1} = find(SessionData.CorrectSide(bTrials) == 1 & sucInd); % correct visual trials
    modIdx{2} = find(SessionData.CorrectSide(bTrials) == 2 & sucInd); % correct audio trials

    animalID{iAnimals} = repmat(iAnimals, 1, size(data.dFOF, 1)); %
    recDepth{iAnimals} = repmat(dataOverview{iAnimals,8}, 1, size(data.dFOF, 1)); %
    recArea{iAnimals} = repmat(dataOverview(iAnimals,7), 1, size(data.dFOF, 1)); %
    
    %make sure trials are evenly sampled from entire session
    binWidth = ceil(length(bTrials) ./ 10);
    edges = 0 : binWidth : binWidth * 10;
            
    [a1, ~, b1] = histcounts(modIdx{1},edges);
    [a2, ~, b2] = histcounts(modIdx{2},edges);
    
    for iBins = 1 : length(a1)
        if a1(iBins) > a2(iBins)               
            cIdx = modIdx{1}(b1 == iBins); %trials in current bin
            b1(ismember(modIdx{1}, cIdx(randperm(length(cIdx), a1(iBins) - a2(iBins))))) = []; %remove trials to make both conditions even in current bin
            modIdx{1}(ismember(modIdx{1}, cIdx(randperm(length(cIdx), a1(iBins) - a2(iBins))))) = []; %remove trials to make both conditions even in current bin
        elseif a1(iBins) < a2(iBins)
            cIdx = modIdx{2}(b2 == iBins); %trials in current bin
            b2(ismember(modIdx{2}, cIdx(randperm(length(cIdx), a2(iBins) - a1(iBins))))) = []; %remove trials to make both conditions even in current bin
            modIdx{2}(ismember(modIdx{2}, cIdx(randperm(length(cIdx), a2(iBins) - a1(iBins))))) = []; %remove trials to make both conditions even in current bin
        end
    end
    fprintf('%d - %d\n', length(modIdx{1}), length(modIdx{2}));
    allData{1,iAnimals} = mean(data.dFOF(:,:,modIdx{1}),3);
    allData{2,iAnimals} = mean(data.dFOF(:,:,modIdx{2}),3);
    
%     cla
%     stdshade((mean(data.dFOF(:,:,modIdx{1}),3)),'r',[],0.5); hold on
%     stdshade((mean(data.dFOF(:,:,modIdx{2}),3)),'g',[],0.5);
%     title([num2str(length(modIdx{1})) ' - ' num2str(length(modIdx{2}))]);
%     pause
end
%%
data1 = cat(1,allData{1,:});
data2 = cat(1,allData{2,:});
cDepth = cat(2,recDepth{:});
cAnimalID = cat(2,animalID{:});
cArea = cat(2,recArea{:});

temp = data1(:, 1 : baseLength) - data2(:, 1 : baseLength);
baseRejIdx = (abs(mean(temp,2)) > 0.05);
rejIdx = ((sum(abs(data1),2) ./ nanstd(sum(abs(data1),2))) > 3) | baseRejIdx; %exclude neurons with weird PSTH

data1(rejIdx,:) = [];
data2(rejIdx,:) = [];
cDepth(rejIdx) = [];
cAnimalID(rejIdx) = [];
cArea(rejIdx) = [];

%%
figure;
areas = {'ALM' 'MM' 'V1' 'RS' 'S1'}; %different recording sites
for iAreas = 1 : length(areas)
    cIdx = ismember(cArea,areas{iAreas});
    
    subplot(1, length(areas), iAreas);   
    
    stdshade(data1(cIdx,1 : baseLength), 'k', (1 : baseLength) / 31, 0.5); hold on
    stdshade(data1(cIdx,baseLength + 2 : size(data1,2)), 'k', (baseLength + 2 : size(data1,2)) / 31, 0.5);
    
    stdshade(data2(cIdx,1 : baseLength), 'r', (1 : baseLength) / 31, 0.5); hold on
    stdshade(data2(cIdx,baseLength + 2 : size(data2,2)), 'r', (baseLength + 2 : size(data2,2)) / 31, 0.5);
    
    title([areas{iAreas} ' - ' num2str(sum(cIdx))]); axis square
    xlim([0 7]);
end




