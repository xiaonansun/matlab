% code to analyze locaNMF results and check changes over time
% cPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\BpodImager\Animals\'; %path to churchlandNAS
cPath = 'Q:\BpodImager\Animals\'; %local path
frameRate = 15;
[dataOverview, ~, ~, ~, segIdx] = rateDiscRecordings;
segFrames = cumsum(floor(segIdx * frameRate)); %max nr of frames per segment
animals = dataOverview(1:10,1);
% animals = {'mSM63' 'mSM64' 'mSM65' 'mSM66' 'Plex01' 'Plex02' 'Fez7' 'Fez10' 'CSP22' 'CSP23'};
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
shrinkMask = arrayResize(allenMask,2) == 1;
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
dimCnt = 20; %number of dimensions from each area that are considered for dimension correlations
groups = {'mSM' 'Fez' 'Plex' 'CSP'};
groupColors = {'r' 'g' 'c' 'k'}; %colors for different groups
animalMarkers = {'+' 'o' '*' 'x'}; %makers for different mice
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
minSize = 1000; %minimum size for areas to use (pixels in both hemispheres)
umapPath = 'R:\BpodImager\umapClust\'; %save path for umap
savePath = 'R:\BpodImager\spatialClust\'; %save path for clustering
cTime = datestr(now,'yyyymmddTHHMM');
areaThresh = 0.75; %threshold to detect locaNMF areas

%% get regions
load([fileparts(cPath(1:end-1)) filesep 'newRegionMap.mat'],'newRegionMap')
[areaMask, areaLabels] = rateDisc_areaMasks(allenMask, minSize); %get different allen areas
areaMask = areaMask(~strcmpi(areaLabels, 'MOB')); %dont use oflactory bulb
areaLabels = areaLabels(~strcmpi(areaLabels, 'MOB'));
sortLabels = {'MOp' 'MOs' 'ACAd' 'PL' 'SSp-n' 'SSp-bfd' 'SSp-ll' ...
              'SSp-m' 'SSp-ul' 'SSp-tr' 'SSp-un' 'SSs' 'VISal' 'VISam' ...
              'VISp' 'VISpm' 'RSPagl' 'RSPd' 'RSPv' 'VISa' 'VISrl'};
              
for iLabels = 1 : length(sortLabels)
	sortIdx(iLabels) = find(ismember(areaLabels, sortLabels{iLabels}));
end
areaMask = areaMask(sortIdx);
areaLabels = areaLabels(sortIdx);

%% run over animals
allA = cell(1,length(animals)); %spatial components
spaceCorr = cell(1,length(animals)); %component spatial correlation maps
tempCorr = cell(1,length(animals)); %component temporal correlation maps
areaCorr = cell(1,length(animals)); %area correlation map
trialCnt = cell(1,length(animals)); %keep trialcount for each recording
iiProb = cell(1,length(animals)); %keep inter-ictal event probability
dimNr = cell(1,length(animals)); %keep number of dimensions in each recording
allAreas = cell(1,length(animals)); %keep area IDs
for iAnimals = 1 : length(animals)
    %current animal
    cAnimal = animals{iAnimals}; % current animal
    recs = rateDisc_getRecsForAnimal(cAnimal, trainingRange, cPath);
    fprintf('Current animal: %s\n', cAnimal);
    
    %go through recordings
    allA{iAnimals} = cell(1, length(recs));
    allAreas{iAnimals} = cell(1, length(recs));
    spaceCorr{iAnimals} = cell(1, length(recs));
    tempCorr{iAnimals} = cell(1, length(recs));
    trialCnt{iAnimals} = NaN(1,length(recs));
    iiProb{iAnimals} = NaN(1,length(recs));
    dimNr{iAnimals} = NaN(1,length(recs));
    areaCorr{iAnimals} = NaN(length(areaMask)*2, length(areaMask)*2, length(recs), 'single'); %sorted correlation map between dimensions
    for iRecs = 1 : length(recs)
        fPath = [cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        try
            load([fPath 'Vc.mat'], 'bTrials');
        catch
            load([fPath 'rsVc.mat'], 'bTrials');
        end
        trialCnt{iAnimals}(iRecs) = length(bTrials);

        load([fPath 'newAC_20_50.mat'], 'A', 'C', 'areas');
        A = A(:,:,areas ~= 1 & areas ~= 255); %don't use OB areas
        C = C(areas ~= 1 & areas ~= 255, :); %don't use OB areas
        areas = areas(areas ~= 1 & areas ~= 255); %don't use OB areas
        A = single(rateDisc_removeOutline(A,10)); %remove some outline to remove variance artifacts
        iiSpikeFrames = findInterictalSpikes(A, C, 2, false); %find interictal spikes
        iiProb{iAnimals}(iRecs) = sum(iiSpikeFrames) / length(iiSpikeFrames); %keep interictal probability
        C = interpOverInterictal(C, iiSpikeFrames); %interpolate over interictal spikes
        
        % compute inter-area correlations
        Cnt = 1;
        cData = NaN(size(C,2), length(areaMask)*2, 'single');
        for x = 1 : length(areaMask)
            cMask = areaMask{x};
            cMask(:, size(cMask,2)/2+1:end) = cMask(:, size(cMask,2)/2+1:end) * 2; %split hemispheres
            cData(:,Cnt) = nanmean(arrayShrink(A,cMask ~= 1, 'merge') * C, 1); %acitivty in current area (left HS)
            cData(:,Cnt+1) = nanmean(arrayShrink(A,cMask ~= 2, 'merge') * C, 1); %acitivty in current area (right HS)
            Cnt = Cnt + 2;
        end
        areaCorr{iAnimals}(:,:,iRecs) = corrcoef(cData); %keep correlation matrix
                
        % compute spatial correlation maps for each component
        A = arrayShrink(arrayResize(A,2), shrinkMask, 'merge');
        covV = cov(C');
        varP1 = dot((A*covV)', A'); % 1 x P
        A(A == 0) = eps(single(1)); %remove absolute zeros
        
        spaceCorr{iAnimals}{iRecs} = NaN(size(A),'single');
        for x = 1 : size(A,2)
            covP = dot((A(:,x)*covV(x,:))', A');
            varP2 = sum((A(:,x) * covV(x,x)) .* A(:,x), 2)'; % 1 x P
            stdPxPy = varP1.^0.5 .* varP2.^0.5; % 1 x P
            spaceCorr{iAnimals}{iRecs}(:,x) = covP./stdPxPy; % 1 x P
        end
        
        % compute temporal correlation maps for each component
        bhvFile = dir([fPath animals{iAnimals} '_SpatialDisc*.mat']);
        load([fPath bhvFile(1).name], 'SessionData');
        bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials
        
        choiceIdx = bhv.ResponseSide == 1; %left response side
        stimIdx = bhv.CorrectSide == 1; %left side is correct
        singleIdx = bhv.DistStim == 0; %detection trials
        
        load([fPath 'opts2.mat'], 'opts');
        opts.frameRate = frameRate;
        load([fPath 'QR.mat'], 'nanIdx');
        Vc = NaN(size(C,1), size(nanIdx,2), 'single');
        Vc(:, ~nanIdx) = C;
        Vc = reshape(Vc, size(C,1), [], length(bTrials));
        Vc = rateDisc_getBhvRealignment(Vc, bhv, segFrames, opts); %align to different trial segments
        
        % get PSTH for each component
%         clear cPSTH
%         cPSTH(:,:,1) = nanmean(Vc,3); %all trial PSTH
%         cPSTH(:,:,2) = nanmean(Vc(:,:,choiceIdx),3); %left choice PSTH
%         cPSTH(:,:,3) = nanmean(Vc(:,:,~choiceIdx),3); %right choice PSTH
%         cPSTH(:,:,4) = nanmean(Vc(:,:,stimIdx & singleIdx),3); %left single stim PSTH
%         cPSTH(:,:,5) = nanmean(Vc(:,:,~stimIdx & singleIdx),3); %right single stim PSTH
        tempCorr{iAnimals}{iRecs} = nanmean(Vc,3)'; %combine PSTHs for clustering
        
        % keep spatial components
        allA{iAnimals}{iRecs} = A;
        allAreas{iAnimals}{iRecs} = areas;
        dimNr{iAnimals}(iRecs) = length(areas);
        
        if any(nanmean(allA{iAnimals}{iRecs},1) == 0)
            disp(recs(iRecs));
            error;
        end
        fprintf('Recording %d/%d\n', iRecs, length(recs));
        clear A C Vc
    end
end
% save([savePath cTime '_areaCorr.mat'],'areaCorr','animals','-v7.3');

%% make area correlation figure
figure('renderer','painters');
cLabels = arrayfun(@(x) [areaLabels{x} '-L'], 1:length(areaLabels), 'UniformOutput', false);
cLabels = [cLabels; arrayfun(@(x) [areaLabels{x} '-R'], 1:length(areaLabels), 'UniformOutput', false)];
cLabels = cLabels(:);
for iGroups = 1 : length(groups)
    
    subplot(2,2,iGroups);
    cGroup = (contains(animals', groups{iGroups}));
    cData = nanmean(cat(3,areaCorr{cGroup}),3);
    ax = imagesc(cData); axis square; 
    caxis([0.5 1])
    
    title(['Predicted region variance - R^2: ' groups{iGroups}]);
    ax.Parent.XTick = 1:size(cData,2);
    ax.Parent.XTickLabelRotation = 45;
    ax.Parent.XTickLabels = cLabels;
    ax.Parent.YTick = 1:size(cData,1);
    ax.Parent.YTickLabels = cLabels;
    ax.Parent.TickLength = [0 0];
    ylabel('Predictor regions');
    niceFigure(ax.Parent)
    axis square; colormap(inferno(256)); colorbar
end

%% correlations across sessions
allCorrs = cat(3,areaCorr{:});
allCorrs = reshape(allCorrs,[],size(allCorrs,3));
allCorrs = corrcoef(allCorrs);
rejIdx = nanmean(allCorrs) < 0.65 | isnan(nanmean(allCorrs));
figure('renderer','painters');
ax = imagesc(allCorrs(~rejIdx,~rejIdx)); caxis([0.5 1]);
colormap(inferno(256)); axis square
[~, b] = cellfun(@size, trialCnt);
aRange = cumsum(b);
for x = 1 : length(aRange)
    aRange(x) = aRange(x) - sum(find(rejIdx) < aRange(x));
end
ax.Parent.TickLength = [0 0];
ax.Parent.YTick = [0 aRange(1:end-1)] + diff([0 aRange]) / 2;
ax.Parent.YTickLabels = animals;
ax.Parent.XTick = [];
nhline(aRange+0.5, 'w', 'linewidth',1)
nvline(aRange+0.5, 'w', 'linewidth',1)
title('Correlation over sessions');
caxis([0.5 0.97])

%% save data for clustering
cLabels = arrayfun(@(x) [areaLabels{x} '-L'], 1:length(areaLabels), 'UniformOutput', false);
cLabels = [cLabels; arrayfun(@(x) [areaLabels{x} '-R'], 1:length(areaLabels), 'UniformOutput', false)];
cLabels = cLabels(:);

allAreaLabels = cat(2,allAreas{:});
allAreaLabels = cat(2,allAreaLabels{:});

animalLabels = zeros(1, sum(cat(2,dimNr{:})), 'single');
groupLabels = zeros(1, sum(cat(2,dimNr{:})), 'single');
Cnt = 0;
for x = 1 : length(dimNr)
    cGroup = find(contains(groups, animals{x}(1:3)));
    animalLabels(Cnt + (1 : sum(dimNr{x}))) = x;
    groupLabels(Cnt + (1 : sum(dimNr{x}))) = cGroup;
    Cnt = Cnt + sum(dimNr{x});
end

% combine components into 2D array
X = cat(2,allA{:});
X = cat(2,X{:});
nanIdx = ~isnan(mean(X,2)); %remove NaN pixels
X = X(nanIdx,:);
nShrinkMask = arrayShrink(nanIdx, shrinkMask, 'split');
nShrinkMask(nShrinkMask == 0) = NaN;
nShrinkMask = isnan(nShrinkMask);
clear allA

spaceX = cat(2,spaceCorr{:});
spaceX = cat(2,spaceX{:});
spaceX = spaceX(nanIdx,:);
clear spaceCorr

tempX = cat(2,tempCorr{:});
tempX = cat(2,tempX{:});
tempX(isnan(mean(tempX,2)),:) = [];

% save data
mkdir([umapPath cTime]);
save([umapPath cTime filesep 'allA.mat'],'X','groupLabels','nShrinkMask','animalLabels','dimNr','allAreaLabels','-v7.3');
save([umapPath cTime filesep 'spaceCorr.mat'],'spaceX','groupLabels','nShrinkMask','animalLabels','dimNr','allAreaLabels','-v7.3');
save([umapPath cTime filesep 'tempCorr.mat'],'tempX','groupLabels','animalLabels','dimNr','allAreaLabels','-v7.3');

%% save data for individual regions
areaMaps = newRegionMap;
areaMaps(:, 1: size(areaMaps)/2) = 256 - areaMaps(:, 1: size(areaMaps)/2);
areaMaps(areaMaps == 256) = 0;

mkdir([umapPath cTime]);
save([umapPath cTime filesep 'areaLabels.mat'],'allAreaLabels','groupLabels','areaMaps','-v7.3');
for iArea = unique(allAreaLabels)
    
    areaIdx = allAreaLabels == iArea;
    areaX = X(:,areaIdx);
    areaX = arrayShrink(areaX, nShrinkMask, 'split');
    
    % get area features like the centroid, orientation etc
    [newX, xArea, xMajor, xMinor, xEccentricity, xOrientation, ...
    xSolidity, xIntensity, xCentroid, xWCentroid] = rateDisc_areaFeatures(areaX, areaThresh);
    
    for x = 1 : sum(areaIdx)
        areaX(:,:,x) = smoothImg(areaX(:,:,x),2,10);
    end
    areaMask = arrayResize(areaMaps == iArea,2) ~= 1 | isnan(mean(areaX,3)); %current region
    areaX = arrayShrink(areaX, areaMask, 'merge');
    areaX = zscore(areaX);
    newX = arrayShrink(newX, areaMask, 'merge');
    binX = newX > 0;
    areaGroups = int16(groupLabels(areaIdx));
    save([umapPath cTime filesep 'areaX_area' num2str(iArea) '.mat'], 'binX', 'newX', 'areaX','areaGroups','nShrinkMask','areaMask', ...
        'xArea','xMajor','xMinor', 'xEccentricity', 'xOrientation', 'xSolidity', 'xIntensity', 'xCentroid', 'xWCentroid','-v7.3');
    
    %spatial correlation maps
    areaX = spaceX(:,areaIdx);
    areaX = arrayShrink(areaX, nShrinkMask, 'split');
    for x = 1 : sum(areaIdx)
        areaX(:,:,x) = smoothImg(areaX(:,:,x));
    end
    areaMask = isnan(mean(areaX,3)); %current region
    areaX = arrayShrink(areaX, areaMask, 'merge');
    save([umapPath cTime filesep 'spaceX_area' num2str(iArea) '.mat'],'areaX','areaGroups','nShrinkMask','areaMask','-v7.3');
end

%% cluster based on spatial features (this is just for visualization)
iArea = 2;
for x = 101 : 500
%         areaIdx = allAreaLabels == iArea;
%     areaX = X(:,areaIdx);
rawFrame = arrayShrink(X(:,x),nShrinkMask,'split');
rawFrame(isnan(rawFrame)) = 0;

% smooth and normalize
cFrame = maxnorm(smoothImg(rawFrame,2,10));

subplot(1,3,1); hold off;
imagesc(maxnorm(rawFrame)); axis image;
hold on; 

% cFrame = imdilate(cFrame,strel('disk',4));
% cFrame = imerode(cFrame,strel('disk',4));
cFrame = imdilate(cFrame > areaThresh,strel('disk',8));
cFrame = imerode(cFrame,strel('disk',8));
cFrame = bwmorph(cFrame,'open'); %break weak connections between different areas
cFrame = bwareaopen(cFrame,50); %don't use areas that are smaller as 10 pixels

areaInfo = regionprops(cFrame, rawFrame, ...
    'Centroid', 'WeightedCentroid', 'MeanIntensity', 'Area', 'MajorAxisLength',  ...
    'MinorAxisLength', 'Eccentricity', 'Orientation', 'Solidity');

if length(areaInfo) > 1
    areas = cat(1,areaInfo(:).Area); %size of different patches
    cIdx = find((areas - (max(areas) / 2)) > 0); %areas that are at least half as large as the biggest one
    
    if length(cIdx) == 1
        areaInfo = areaInfo(cIdx);
    else
        areaInfo = [];
    end
end

contour(cFrame,'w');

subplot(1,3,2);
imagesc(cFrame); axis image;
title(sum(cFrame(:) > areaThresh));

rawFrame(~cFrame) = 0;
subplot(1,3,3);
imagesc(maxnorm(rawFrame)); axis image; colormap(viridis(256));
if ~isempty(areaInfo)
    plotEllipse(gca, areaInfo); hold off;
end

pause
end

%% show umap result for spatial components and color in regions
% clustFile = '20200612T1221_umap_group_1_area_254.mat';
% load([umapPath clustFile],'spatialUmap')
% newAreaLabels = removeGaps(allAreaLabels); %adjust area indicies
% 
% figure
% subplot(3,2,[1 3 5])
% ax = scatter(spatialUmap(:,1), spatialUmap(:,2), 25, newAreaLabels(groupLabels == 1), '.', 'Linewidth',2); axis square
% ax.Parent.Visible = 'off';
% title('umap - all Areas');
% 
% % show examples of invidiual clusters
% subplot(3,2,2); x = 2;
% cFile = strrep(clustFile,'allA',['group_1_area_' num2str(x) '.mat']);
% load([umapPath cFile], 'spatialLabels')
% cIdx = find(allAreaLabels == x & groupLabels == 1);
% cMap = zeros(sum(~shrinkMask(:)),1);
% cMap(nanIdx) = nanmean(X(:,cIdx(spatialLabels == 6)),2);
% imageScale(arrayShrink(cMap,shrinkMask,'split')); 
% rateDisc_plotAllenOutline(gca); caxis([0 0.75]);
% colormap(inferno);
% title(['example - area ' num2str(x)]);
% 
% subplot(3,2,4); x = 4;
% cFile = strrep(clustFile,'allA',['group_1_area_' num2str(x) '.mat']);
% load([umapPath cFile], 'spatialLabels')
% cIdx = find(allAreaLabels == x & groupLabels == 1);
% cMap = zeros(sum(~shrinkMask(:)),1);
% cMap(nanIdx) = nanmean(X(:,cIdx(spatialLabels == 2)),2);
% imageScale(arrayShrink(cMap,shrinkMask,'split')); 
% rateDisc_plotAllenOutline(gca); caxis([0 0.9]);
% colormap(inferno);
% title(['example - area ' num2str(x)]);
% 
% subplot(3,2,6); x = 8;
% cFile = strrep(clustFile,'allA',['group_1_area_' num2str(x) '.mat']);
% load([umapPath cFile], 'spatialLabels')
% cIdx = find(allAreaLabels == x & groupLabels == 1);
% cMap = zeros(sum(~shrinkMask(:)),1);
% cMap(nanIdx) = nanmean(X(:,cIdx(spatialLabels == 2)),2);
% imageScale(arrayShrink(cMap,shrinkMask,'split')); 
% rateDisc_plotAllenOutline(gca); caxis([0 0.9]);
% title(['example - area ' num2str(x)]);
% colormap(inferno);

%% check umap clusters and merge if needed
% load results and take a look
for iArea = unique(allAreaLabels)
    
    % get clustering
    fIdx = find(allAreaLabels == iArea);
    cFile = [umapPath cTime filesep 'umap_area_' num2str(iArea) '.mat'];
    load(cFile,'binUmap','binLabels','cIdx')
    
    cFile = [umapPath cTime filesep 'areaX_area' num2str(iArea) '.mat'];
    load(cFile,'binX','areaMask')
    binX = binX(:,cIdx);

    temp = find(cIdx);
    temp(binLabels > -1) = [];
    cIdx(temp) = false;
    fIdx(~cIdx) = [];
    binUmap(binLabels == -1,:) = [];
    binX(:,binLabels == -1,:) = [];
    binLabels(binLabels == -1) = [];
    
    figure
    subplot(1,2,1);
    scatter(binUmap(:,1),binUmap(:,2),25,binLabels,'o','LineWidth',5); axis square
    
    % check cluster
    outFile = [umapPath cTime filesep 'umap_area_' num2str(iArea) '_checked.mat'];

    % check best cluster for each component
    a = corrcoef(X(:,fIdx));
    cLabels = unique(spatialLabels);
    tempLabels = spatialLabels;
    ownCorr = NaN(1,length(spatialLabels));
    for x = 1 : length(spatialLabels)
        ownCorr(x) = nanmean(a(x,spatialLabels == spatialLabels(x))); %correlation to own cluster
        for y = 1:length(cLabels)
            if nanmean(a(x,spatialLabels == cLabels(y))) > ownCorr(x) %more correlated to different cluster
                ownCorr(x) = nanmean(a(x,spatialLabels == cLabels(y)));
                fprintf('Switch from %d to %d\n',spatialLabels(x),cLabels(y));
                tempLabels(x) = cLabels(y);
            end
        end
    end
    spatialLabels = tempLabels;
    
    %remove too weak cluster members
    spatialLabels(ownCorr < 0.8) = [];
    spatialUmap(ownCorr < 0.8,:) = [];
    cIdx(ownCorr < 0.8) = [];
    
    subplot(1,2,2);
    scatter(spatialUmap(:,1),spatialUmap(:,2),25,spatialLabels','o','LineWidth',5); axis square
    
    %manually check all clusters and merge similar clusters together
    spatialLabels = removeGaps(spatialLabels); %adjust cluster indicies
%     rateDisc_checkCluster(spatialLabels', X(:,cIdx)',nShrinkMask,4,outFile);
    rateDisc_checkCluster(removeGaps(binLabels)', rX',nShrinkMask,4,outFile);
    
    %%
    figure
    subplot(1,2,1);
    scatter(spatialUmap(:,1),spatialUmap(:,2),25, spatialLabels,'o'); axis square
    subplot(1,2,2);
    scatter(spatialUmap(spatialLabels>-1,1),spatialUmap(spatialLabels>-1,2),25,spatialLabels(spatialLabels>-1),'o'); axis square
    
    cFile2 = [umapPath 'checked' filesep strrep(cFile1, '.mat', '_checked.mat')];
    rateDisc_checkCluster(spatialLabels'+2, X(:,cIdx)',nShrinkMask,4,cFile2);

end

% check results and take a look
checkFile = dir([umapPath '*checked.mat']);
checkFile = {checkFile.name};
for iArea = 254
% for iArea = unique(areaLabels)
    
    % get clustering    
    cFile = [umapPath clustFile{contains(clustFile, ['area' int2str(iArea) '_clustered.mat'])}]; %get file
    load(cFile); disp(cFile);
    cFile = [umapPath checkFile{contains(checkFile, ['area' int2str(iArea) '_checked.mat'])}]; %get file
    load(cFile); disp(cFile);
    
    figure
    scatter(umapOut(hdbLabels>-1,1),umapOut(hdbLabels>-1,2),25,hdbLabels(hdbLabels>-1),'o'); axis square
    scatter(umapOut(hdbLabels>-1,1),umapOut(hdbLabels>-1,2),25,newT(hdbLabels>-1),'o'); axis square

    cFile = strrep(cFile, 'clustered.mat', 'checked.mat');
    rateDisc_checkCluster(newT,X(areaLabels == iArea,:),nShrinkMask,4,cFile);

end
% areaY = tsne(areaX,'Algorithm','barneshut','Perplexity',pplx);

%%
figure; title(['area: ' num2str(iArea) ' - Peplexity: ' num2str(pplx)]); hold on
for x = 1 : length(animals)
    cGroup = contains(groups, animals{x}(1:3));
    cAnimal =  contains(animals, groups{cGroup});
    cAnimal = contains(animals(cAnimal)',animals{x});
    cIdx = animalLabels(areaIdx) == x; 
    plot(R(cIdx,2),R(cIdx,1), 'color', groupColors{cGroup}, 'Marker',animalMarkers{cAnimal}, 'Linestyle', 'none','LineWidth',2)
    axis square; 
end

% %%
% for x = 1:4
%     figure
%     for y = 1 : 4
%         subplot(2,2,y)
%         cIdx = ismember(groupLabels,[x y]); gscatter(Y(cIdx,1),Y(cIdx,2),groupLabels(cIdx)'); axis square; legend('off')
%     end
% end

%% spatial clustering
Y1 = round((Y - min(Y))*10)+1;
cFrame = false(max(Y1));
xx = sub2ind(size(cFrame),Y1(:,1),Y1(:,2));

cFrame(xx) = true;
cFrame = imdilate(cFrame,strel('disk',7));
cFrame = imerode(cFrame,strel('disk',4));
cFrame = imclose(cFrame,strel('disk',4));
cFrame = bwareaopen(cFrame,10); %don't use areas that are smaller as 10 pixels
imagesc(cFrame); axis image;

figure;
imagesc(bwlabel(cFrame)); axis image; hold on
areaInfo = regionprops(cFrame, 'Area', 'Centroid', 'MajorAxisLength','MinorAxisLength','Orientation','Solidity');
cFrame = bwlabel(cFrame);
clustID = cFrame(xx);

for iRuns = 1 : length(areaInfo)
    phi = linspace(0,2*pi,50);
    cosphi = cos(phi);
    sinphi = sin(phi);
    
    theta = pi*areaInfo(iRuns).Orientation/180;
    R = [ cos(theta)   sin(theta)
        -sin(theta)   cos(theta)];
    
    xy = [(areaInfo(iRuns).MajorAxisLength/2)*cosphi; (areaInfo(iRuns).MinorAxisLength/2)*sinphi];
    xy = R*xy;
    x = xy(1,:) + areaInfo(iRuns).Centroid(1);
    y = xy(2,:) + areaInfo(iRuns).Centroid(2);
    
    plot(x,y,'r','LineWidth',2);
    plot(areaInfo(iRuns).Centroid(1),areaInfo(iRuns).Centroid(2),'xr','LineWidth',2);
    plot(Y1(clustID == iRuns ,2),Y1(clustID == iRuns,1),'o','LineWidth',2)
end


%% kmeans overclustering with manual merging step
areas = unique(areaLabels);
for iAreas = areas
areaIdx = areaLabels == iAreas;
nrRows = 4;
T = kmeans(X(areaIdx,:),20);
cFile = [savePath cTime '_area' num2str(iAreas) '_kmeans.mat'];
rateDisc_checkCluster(newT,X(areaIdx,:),nShrinkMask,nrRows,cFile);
end
