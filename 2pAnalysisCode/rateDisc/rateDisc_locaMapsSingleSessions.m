% code to analyze locaNMF results and check changes over time
% cPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\BpodImager\Animals\'; %path to churchlandNAS
cPath = 'Q:\BpodImager\Animals\'; %local path
dataOverview = rateDiscRecordings;
animals = dataOverview(1:10,1);
% animals = {'mSM63' 'mSM64' 'mSM65' 'mSM66' 'Plex01' 'Plex02' 'Fez7' 'Fez10' 'CSP22' 'CSP23'};
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
shrinkMask = arrayResize(allenMask,2) == 1;
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
dimCnt = 20; %number of dimensions from each area that are considered for dimension correlations
maxCluster = 100; %how many spatial clusters should be used
groups = {'mSM' 'Fez' 'Plex' 'CSP'};
groupColors = {'r' 'g' 'c' 'k'}; %colors for different groups
animalMarkers = {'+' 'o' '*' 'x'}; %makers for different mice
savePath = 'Q:\BpodImager\spatialClust\'; %save path for clustering
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS

%% get regions
load([fileparts(cPath(1:end-1)) filesep 'regionMap.mat'],'regionMap')
load([fileparts(cPath(1:end-1)) filesep 'largeRegionMap.mat'],'largeRegionMap')
regLabels = {'OB' 'M2' 'M1' 'SSb' 'SSf' 'Aud' 'PPC' 'RS' 'V2' 'V1'}; %labels for different regions
regs = 1:max(regionMap(:)); %regions based on 'regionMap'
% regs = 1:max(largeRegionMap(:)); %regions based on 'regionMap'
regLabels = regLabels(1:length(regs));

%% run over animals
dimCorr = cell(1,length(animals)); %sorted correlation map between components
regCorr = cell(1,length(animals)); %sorted correlation map between areas
dimNr = cell(1,length(animals));
iiProb = cell(1,length(animals));
trialCnt = cell(1,length(animals));
varMaps = cell(1,length(animals));
allA = cell(1,length(animals));

for iAnimals = 1 : length(animals)
    %current animal
    cAnimal = animals{iAnimals}; % current animal
    recs = rateDisc_getRecsForAnimal(cAnimal, trainingRange, cPath);
    fprintf('Current animal: %s\n', cAnimal);
    
    %pre-allocate
    dimCorr{iAnimals} = NaN(length(regs) * dimCnt*2, length(regs) * dimCnt*2, length(recs), 'single'); %sorted correlation map between dimensions
    regCorr{iAnimals} = NaN(length(regs), length(regs), length(recs), 'single'); %sorted correlation map between dimensions
    varMaps{iAnimals} = NaN(sum(~allenMask(:)), length(recs), 'single'); % variance maps
    allA{iAnimals} = cell(1,length(recs)); % variance maps
    dimNr{iAnimals} = NaN(1, length(recs)); %get number of components for each recording
    iiProb{iAnimals} = NaN(1, length(recs)); %get probability of interictal events in each recording
    trialCnt{iAnimals} = NaN(1, length(recs)); %get number of trials in each recording
    %go through recordings
    for iRecs = 1 : length(recs)
        fPath = [cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        try
            load([fPath 'Vc.mat'], 'bTrials');
        catch
            load([fPath 'rsVc.mat'], 'bTrials');
        end
        trialCnt{iAnimals}(iRecs) = length(bTrials);

        load([fPath 'AC.mat'], 'A', 'C', 'areas');
%         load([fPath 'largeAC.mat'], 'A', 'C', 'areas');
%         load([fPath 'largeAC_20_50.mat'], 'A', 'C', 'areas');
        A = single(rateDisc_removeOutline(A,3)); %remove some outline to remove variance artifacts

        iiSpikeFrames = findInterictalSpikes(A, C, 2, false); %find interictal spikes
        iiProb{iAnimals}(iRecs) = sum(iiSpikeFrames) / length(iiSpikeFrames); %keep interictal probability
        C = interpOverInterictal(C, iiSpikeFrames); %interpolate over interictal spikes
        
        % compute total variance map
        [~, varP1] = rateDisc_modelCorr(C, C, A); 
        varP1 = arrayShrink(varP1,isnan(A(:,:,1)), 'split'); 
        varMaps{iAnimals}(:,iRecs) = arrayShrink(varP1,allenMask, 'merge'); 
        
        % compute inter-region correlations
        [regCorr{iAnimals}(:,:,iRecs), dimCorr{iAnimals}(:,:,iRecs)] = rateDisc_regionCorrLocal(A, C, areas, regs, dimCnt);
        
        % downsample and select components to keep
        A = A(:,:,areas ~= 1 & areas ~= 255); %don't use OB areas
        allA{iAnimals}{iRecs} = arrayShrink(arrayResize(A,2), shrinkMask, 'merge');
        areas = areas(areas ~= 1 & areas ~= 255); %don't use OB areas
        allAreas{iAnimals}{iRecs} = areas;
        dimNr{iAnimals}(iRecs) = size(A,3);
        
        if any(nanmean(allA{iAnimals}{iRecs},1) == 0)
            disp(recs(iRecs));
            error;
        end
        fprintf('Recording %d/%d\n', iRecs, length(recs));
        clear A
    end
end
save([savePath datestr(now,'yyyymmddTHHMM') '_varMaps.mat'],'varMaps','animals','-v7.3');
save([savePath datestr(now,'yyyymmddTHHMM') '_regCorr.mat'],'regCorr','regLabels','animals','groups','-v7.3');

%% show variance maps
figure('renderer','painter');
cRanges = [2.5 1.5 5.5 6]*10^-4;
gPos = cell(1,length(groups));
for iGroups = 1 : length(groups)
    
    Cnt = 0;
    cGroup = find(contains(animals, groups{iGroups}))';
    gPos{iGroups} = cGroup + iGroups;
    for x = 1 : length(cGroup)
        
        Cnt = Cnt+1;
        cIdx = zscore(nanmean(varMaps{cGroup(x)},1)) < 3.5;
        meanVar{cGroup(x)} = nanmean(varMaps{cGroup(x)}(:,cIdx),1);
        groupId{cGroup(x)} = ones(1,sum(cIdx)).*(cGroup(x) + iGroups);
        cData = arrayShrink(varMaps{cGroup(x)}, allenMask,'split');
        cData = nanmean(rateDisc_removeOutline(cData,3),3); %remove some outline to remove variance artifacts
        
        subplot(length(groups),5,(iGroups-1)*5 + Cnt);
        mapImg = imageScale(cData); colormap('inferno'); drawnow; title([groups{iGroups} '-' int2str(x)]);
        
        hold(mapImg.Parent, 'on'); caxis([0 cRanges(iGroups)])
        for y = 1: length(dorsalMaps.edgeOutlineSplit)
            plot(dorsalMaps.edgeOutlineSplit{y}(:,2), dorsalMaps.edgeOutlineSplit{y}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
    end
    
    cData = cat(2,varMaps{cGroup});
    cData = arrayShrink(nanmean(cData,2), allenMask,'split');
    
    subplot(length(groups),5,(iGroups-1)*5 + Cnt + 1);
    mapImg = imageScale(cData); colormap('inferno'); drawnow; title([groups{iGroups} '-Mean']);
    
    hold(mapImg.Parent, 'on'); caxis([0 cRanges(iGroups)])
    for y = 1: length(dorsalMaps.edgeOutlineSplit)
        plot(dorsalMaps.edgeOutlineSplit{y}(:,2), dorsalMaps.edgeOutlineSplit{y}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% show boxplot   
figure('renderer','painters');
boxplot(cat(2,meanVar{:}), cat(2,groupId{:}), 'positions', cat(2,gPos{:}), 'labels', animals)
title('cross-session variance per animal');
ylabel('variance'); niceFigure(gca); axis square;
ax = gca;
ax.XTickLabelRotation = 45;
ax.YLim(2) = sum(ax.YLim);
ax.YLim(1) = 0;

figure('renderer','painters');
hold on;
for x = 1 : length(animals)
    plot(smooth(meanVar{x}));
end
niceFigure(gca); axis square;

%% make regional prediction figure
clear sortIdx
sortLabels = {'M2' 'M1' 'SSb' 'SSf' 'Aud' 'RS' 'PPC' 'V2' 'V1'}; %labels for different regions - use this to change order of regions
for iLabels = 1 : length(sortLabels)
    sortIdx(iLabels) = find(ismember(regLabels, sortLabels{iLabels}));
end

figure;
for iGroups = 1 : length(groups)
    subplot(2,2,iGroups);
    cGroup = (contains(animals', groups{iGroups}));
    cData = nanmean(cat(3,regCorr{cGroup}),3);
    ax = imagesc(cData(sortIdx,sortIdx)); axis square; 
    caxis([0 0.6])
    
    title(['Predicted region variance - R^2: ' groups{iGroups}]);
    ax.Parent.XTick = 1:size(cData,2);
    ax.Parent.XTickLabelRotation = 45;
    ax.Parent.XTickLabels = sortLabels;
    ax.Parent.YTick = 1:size(cData,1);
    ax.Parent.YTickLabels = sortLabels;
    ax.Parent.TickLength = [0 0];
    ylabel('Predictor regions');
    niceFigure(ax.Parent)
    axis square; colormap(parula(256)); colorbar
end
   
%% do clustering stuff
areaLabels = cat(2,allAreas{:});
areaLabels = cat(2,areaLabels{:});

animalLabels = zeros(1, sum(cat(2,dimNr{:})), 'single');
groupLabels = zeros(1, sum(cat(2,dimNr{:})), 'single');
Cnt = 0;
for x = 1 : length(dimNr)
    cGroup = find(contains(groups, animals{x}(1:3)));
    animalLabels(Cnt + (1 : sum(dimNr{x}))) = x;
    groupLabels(Cnt + (1 : sum(dimNr{x}))) = cGroup;
    Cnt = Cnt + sum(dimNr{x});
end

X = cat(2,allA{:});
X = cat(2,X{:});

nanIdx = ~isnan(mean(X,2)); %remove NaN pixels
X = X(nanIdx,:)';
nShrinkMask = arrayShrink(nanIdx, shrinkMask, 'split');
nShrinkMask(nShrinkMask == 0) = NaN;
nShrinkMask = isnan(nShrinkMask);
% save([savePath datestr(now,'yyyymmddTHHMM') '_largeAC_20_50.mat'],'X','-v7.3');
% save([savePath datestr(now,'yyyymmddTHHMM') '_otherAC_20_50.mat'],'areas','areaLabels','nShrinkMask','animalLabels','groupLabels','dimNr','pplx','-v7.3');
save([umapPath datestr(now,'yyyymmddTHHMM') '_umap_all.mat'],'X','groupLabels','nShrinkMask','animalLabels','groupLabels','dimNr','areaLabels','-v7.3');

%% tsne
pplx = 200;
Y = tsne(X,'Algorithm','barneshut','Perplexity',pplx);
save([savePath datestr(now,'yyyymmddTHHMM') '_largeAC_20_50_tsne_all.mat'],'Y','areas','areaLabels','nShrinkMask','animalLabels','groupLabels','dimNr','pplx','-v7.3');

%% show result
figure; title(['Peplexity: ' num2str(pplx)]); hold on
for x = 1 : length(animals)
    
    cGroup = contains(groups, animals{x}(1:3));
    cAnimal =  contains(animals, groups{cGroup});
    cAnimal = contains(animals(cAnimal)',animals{x});
    cIdx = animalLabels == x; 
    plot(Y(cIdx,2),Y(cIdx,1), 'color', groupColors{cGroup}, 'Marker',animalMarkers{cAnimal}, 'Linestyle', 'none','LineWidth',2)
    axis square; 
    
end

%%umap clustering
umapPath = 'Q:\BpodImager\umapClust\';
save([umapPath 'allLabels.mat'], 'areaLabels', 'animalLabels', 'groupLabels')
lHS = nShrinkMask;
lHS(:, round(size(nShrinkMask,2)/2):end) = true;
rHS = nShrinkMask;
rHS(:, 1 : round(size(nShrinkMask,2)/2)) = true;

% save areas
% for iArea = unique(areaLabels)
%     
%     areaIdx = areaLabels == iArea;
%     areaX = X(areaIdx, :);
%     areaX = arrayShrink(areaX',nShrinkMask,'split');
%     if iArea > 2^7
%         areaX = arrayShrink(areaX,lHS,'merge');
%     else
%         areaX = arrayShrink(areaX,rHS,'merge');
%     end
%     
%     areaGroups = int16(groupLabels(areaIdx));
%     save([umapPath datestr(now,'yyyymmddTHHMM') '_umap_area' num2str(iArea) '.mat'],'areaX','areaGroups','nShrinkMask','lHS','rHS','-v7.3');
% end

% load results and take a look
clustFile = dir([umapPath '*clustered.mat']);
clustFile = {clustFile.name};
for iArea = 2
% for iArea = unique(areaLabels)
    
    % get clustering
    cFile = [umapPath clustFile{contains(clustFile, ['area' int2str(iArea) '_clustered.mat'])}]; %get file
    load(cFile); disp(cFile);
    
    figure
    scatter(umapOut(hdbLabels>-1,1),umapOut(hdbLabels>-1,2),25,hdbLabels(hdbLabels>-1),'o'); axis square
    
    cFile = strrep(cFile, 'clustered.mat', 'checked.mat');
    rateDisc_checkCluster(hdbLabels'+2,X(areaLabels == iArea,:),nShrinkMask,4,cFile);

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
cFile = [savePath datestr(now,'yyyymmddTHHMM') '_area' num2str(iAreas) '_kmeans.mat'];
rateDisc_checkCluster(newT,X(areaIdx,:),nShrinkMask,nrRows,cFile);
end

%% check clusters
% figure
% for iRuns = 1 : max(clustID)
%     cIdx = groupLabels' == 1 & clustID == iRuns;
%     subplot(1,2,1);
%     imagesc(corrcoef(X(cIdx,:)')); axis image
%     subplot(1,2,2);
%     imagesc(arrayShrink(nanmean(X(cIdx,:),1)',nShrinkMask,'split')); axis image
%     pause;
% end

%% hirarchical clustering
% smallY = tsne(X(areaLabels == 2,:),'Algorithm','barneshut');

% Y = pdist(X(areaLabels == 2,:),'euclidean');
% Z = linkage(Y,'ward');
% leafOrder = optimalleaforder(Z,Y);

% T = cluster(Z,'Depth',maxCluster,'CutOff',5);
% T = cluster(Z,'MaxClust',24);
% dendrogram(Z,length(unique(T)),'Reorder',leafOrder)
% I = inconsistent(Z,maxCluster); %4th column is the inconsistent value



%% compare correlation matrices
% cData = cat(3,regCorr{:});
% cData = cData(sortIdx,sortIdx,:);
% cData = arrayShrink(cData,diag(true(1,length(sortIdx))), 'merge');
% a = corrcoef(cData).^2;
% b = tril(a,-1);
% 
% figure
% imagesc(a);
% for x = 1 : length(animals)
%     recCnt(x) = size(regCorr{x},3);
% end
%     
% sameCorr = cell(1,length(animals));
% for iAnimals = 1 : length(animals)
%     if iAnimals == 1
%         cIdx = 1 : recCnt(1);
%     else
%         cIdx = sum(recCnt(1:iAnimals-1))+1 : sum(recCnt(1:iAnimals));
%     end
%     c = b(cIdx,cIdx);
%     sameCorr{iAnimals} = c(c ~= 0);
% end
% 
% 
% Cnt = 0;
% for iGroups = 1 : length(groups)
% 
%     cGroup = find(contains(animals', groups{iGroups}));
%     recCnt = NaN(1, length(cGroup));
%     for x = 1 : length(cGroup)
%         recCnt(x) = size(regCorr{x},3);
%     end
%     
%     % compute variance across animals of the same group
%     cData = cat(3,regCorr{cGroup});
%     cData = cData(sortIdx,sortIdx,:);
%     cData = arrayShrink(cData,diag(true(1,length(sortIdx))), 'merge');
%     a = corrcoef(cData).^2;
%     b = tril(a,-1);
%     
%     % get variance for each animal
%     for x = 1 : length(cGroup)
%         if x == 1
%             cIdx = 1 : recCnt(1);
%         else
%             cIdx = sum(recCnt(1:x-1))+1 : sum(recCnt(1:x));
%         end
%         
%         %same animal, across sessions
%         Cnt = Cnt + 1;
%         c = b(cIdx,cIdx);
%         sameCorr{Cnt} = c(c ~= 0);
% 
%         %same animal, across sessions
% 
%         
%     end
% end

%%
% figure
% for iGroups = 1 : length(groups)
%     cGroup = (contains(animals', groups{iGroups}));
%     subplot(2,2,iGroups);
%     ax = imagesc(nanmean(real(regCorr(:,:,cGroup)),3)); 
%     title(['Predicted variance - R^2: ' groups{iGroups}]);
%     ax.Parent.XTick = 1:size(regCorr,2);
%     ax.Parent.XTickLabelRotation = 45;
%     ax.Parent.YTickLabels = regLabels;
%     ax.Parent.YTick = 1:size(regCorr,2);
%     ax.Parent.TickLength = [0 0];
%     ylabel('Predictor regions');
%     niceFigure(ax.Parent)
%     axis square; colormap('inferno');
%     caxis([0.1 0.9]);
%     hline([2.5 6.5], '--w') %outline for motor/ss regions
%     vline([2.5 6.5], '--w') %outline for motor/ss regions
%     hline([7.5 10.5], '--g') %outline for vision regions
%     vline([7.5 10.5], '--g') %outline for vision regions
%     niceFigure(ax.Parent)
%     if iGroups > 2
%         ax.Parent.XTickLabels = regLabels;
%         xlabel('Predicted regions');
%     else
%         ax.Parent.XTickLabels = [];
%     end
% end

%% make component correlation maps
% lineIdx = repmat(allDims',dimCnt,1);
% lineIdx = reshape(lineIdx, 1, []);
% a = mean(dimCorr,3);
% rejIdx = isnan(nanmean(a,2));
% lineIdx(rejIdx) = [];
% a(rejIdx,:) = [];
% a(:,rejIdx) = [];
% 
% clear regMean regOff
% for iRegs = 1 : length(reg)
%     if any(ismember(lineIdx, reg{iRegs}))
%         regMean(iRegs) = nanmean(find(ismember(lineIdx, reg{iRegs})));
%         regOff(iRegs) = find(ismember(lineIdx, reg{iRegs}),1, 'last');
%     else
%         regMean(iRegs) = NaN;
%     end
% end
% 
% figure
% b = tril(smoothImg(nanmean(real((a)),3),2,1.1));
% b(b == 0) = NaN;
% ax = imagesc(b); title('Predicted variance - R^2');
% ax.Parent.XTick = regMean;
% ax.Parent.XTickLabelRotation = 45;
% ax.Parent.XTickLabels = regLabels(~isnan(regMean));
% ax.Parent.YTickLabels = regLabels(~isnan(regMean));
% ax.Parent.YTick = regMean;
% ax.Parent.TickLength = [0 0];
% set(ax,'AlphaData',~isnan(ax.CData)); %make NaNs transparent.
% xlabel('Predicted regions');
% ylabel('Predictor regions');
% axis square; colormap('parula');
% caxis([-0.03 0.03]);
% hline(regOff(1:end-1)+0.5, '--w');
% vline(regOff(1:end-1)+0.5, '--w');
% niceFigure(ax.Parent)
% 
% %% check similarity of spatial components
% cGroup = find(contains(animals', groups{1}));
% a = cat(2,allA{:});
% b = a(~isnan(mean(a,2)),:);
% b = corrcoef(b(:,nanmean(b,1) ~= 0));
% 
% %%
% figure;
% ax = imagesc(b(:,1:sum(dimNr(1)))'); colorbar; axis image;
% % ax = imagesc(b(:,sum(dimNr(1)):sum(dimNr(1:2)))'); colorbar; axis image;
% title('Correlation across animals');
% ylabel(['Spatial components - ' animals{cGroup(1)}]);
% xlabel(['Spatial components - ' groups{1}]);
% niceFigure(ax.Parent);
