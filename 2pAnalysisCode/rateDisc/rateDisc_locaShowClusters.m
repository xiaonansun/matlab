% function rateDisc_locaShowClusters(cPath)

load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
dataOverview = rateDiscRecordings;
groups = {'mSM' 'Fez' 'Plex' 'CSP'};
animals = dataOverview(1:10,1);
clustPath = 'R:\BpodImager\spatialClust\';
umapPath = 'R:\BpodImager\umapClust\';
clustThresh = 1/2; %threshold for cluster specificity

% load data
load([umapPath '20200512T0009_umap_all.mat'],'X','dimNr','nShrinkMask','shrinkMask');
load([clustPath '20200511T2312_varMaps.mat'],'varMaps');
load([clustPath '20200511T2312_regCorr.mat'],'regCorr','regLabels','animals','groups');

load('Q:\BpodImager\Animals\mSM66\SpatialDisc\11-Jun-2018\newAC_20_50.mat','A');
A = arrayShrink(A,allenMask,'merge');
A = arrayShrink(A,allenMask,'split');

%% show locanMF results
load([fileparts(cPath(1:end-1)) filesep 'newRegionMap.mat'],'newRegionMap')
regionMap = newRegionMap;
lHS = ~allenMask;
lHS(:, 1: round(size(allenMask,2)/2)+1) = false;
regionMapPlot = regionMap;
regionMapPlot = arrayShrink(regionMapPlot,allenMask,'merge');
regionMapPlot = arrayShrink(regionMapPlot,allenMask,'split');
regionMapPlot(regionMapPlot == 0) = NaN;
regionMapPlot(lHS) = regionMapPlot(lHS)+max(regionMapPlot(:));
for iRegions = (1:10) + max(regionMap(:))
    cFrame = regionMapPlot == iRegions;
    cFrame = imclose(cFrame,strel('disk',6));
    regionMapPlot(cFrame) = iRegions;
end

figure('renderer','painters');
subplot(1,4,1);
cImg = imageScale(regionMapPlot);
rateDisc_plotAllenOutline(gca,'L');
colormap(cImg.Parent,hsv(256));
caxis([0 max(regionMapPlot(:))]); hold on;
for iRegions = (1:10) + max(regionMap(:))
    cFrame = regionMapPlot == iRegions;
    areaOutline= bwboundaries(cFrame); %outline of selected area
    for x = 1 : length(areaOutline)
        plot(areaOutline{x}(:,2),areaOutline{x}(:,1),'k')
    end
end

subplot(1,4,2);
cImg = imageScale(A(:,:,51));
caxis([0 0.9]); 
colormap(cImg.Parent,inferno(256));
rateDisc_plotAllenOutline;

subplot(1,4,3);
cImg = imageScale(A(:,:,40));
caxis([0 0.9]); 
colormap(cImg.Parent,inferno(256));
rateDisc_plotAllenOutline;

subplot(1,4,4);
cImg = imageScale(A(:,:,43));
caxis([0 0.9]); 
colormap(cImg.Parent,inferno(256));
rateDisc_plotAllenOutline;


%% show umap result
% get all clusters
umapFile = dir([umapPath '*umap_all_clustered*']);
load([umapPath umapFile.name], 'umapOut', 'hdbLabels')
load([umapPath 'allLabels.mat'], 'allAreaLabels', 'animalLabels', 'groupLabels')

figure('renderer','painters'); 
hold on; title('umap, all PCs');
for x = 1 : length(animals)
    cGroup = contains(groups, animals{x}(1:3));
    cAnimal = contains(animals, groups{cGroup});
    cAnimal = contains(animals(cAnimal)',animals{x});
    cIdx = animalLabels == x; 
    plot(umapOut(cIdx,2),umapOut(cIdx,1), 'color', groupColors{cGroup}, 'Marker','.', 'MarkerSize',10, 'Linestyle', 'none','LineWidth',2)
%     plot(umapOut(cIdx,2),umapOut(cIdx,1), 'color', groupColors{cGroup}, 'Marker',animalMarkers{cAnimal}, 'Linestyle', 'none','LineWidth',2)
    axis square; 
end

%% make regional prediction figure
clear sortIdx
sortLabels = {'M2' 'M1' 'SSb' 'SSf' 'Aud' 'RS' 'PPC' 'V2' 'V1'}; %labels for different regions - use this to change order of regions
for iLabels = 1 : length(sortLabels)
    sortIdx(iLabels) = find(ismember(regLabels, sortLabels{iLabels}));
end

% functional correlations across regions for each group
figure('rendere','painters');
for iGroups = 1 : length(groups)
    subplot(2,2,iGroups);
    cGroup = (contains(animals', groups{iGroups}));
    cData = nanmean(cat(3,regCorr{cGroup}),3);
    ax = imagesc(cData(sortIdx,sortIdx));
    caxis([0 0.6]);
    
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

% functional correlations across regions for each group
% groupOrder = [1:4 7:8 5:6 9:10];
groupOrder = 1:10;
cData = cat(3,regCorr{groupOrder});
cData = cData(sortIdx,sortIdx,:);
cData = arrayShrink(cData,diag(true(1,length(sortIdx))), 'merge');
a = corrcoef(cData);
% b = tril(a,-1);

figure('renderer','painters');
imagesc(smoothImg(a,2,1));
% Cnt = 1;
% for x = 1 : length(groups)
%     cGroup = (contains(animals(groupOrder)', groups{x}));
%     nhline(Cnt + length(cat(2,dimNr{cGroup}))-0.5,'w','linewidth',2);
%     nvline(Cnt + length(cat(2,dimNr{cGroup}))-0.5,'w','linewidth',2);
%     Cnt = Cnt + length(cat(2,dimNr{cGroup}));
% end
Cnt = 1;
for x = 1 : length(animals)
    nhline(Cnt + length(cat(2,dimNr{groupOrder(x)})) - 0.5,'w','linewidth',2);
    nvline(Cnt + length(cat(2,dimNr{groupOrder(x)})) - 0.5,'w','linewidth',2);
    Cnt = Cnt + length(cat(2,dimNr{groupOrder(x)}));
end
axis square; colormap viridis;
caxis([0.2 .8]);
 
%% check number of clusters for each recording
plotUMAP = false; %flag to show umap clusters
load([umapPath 'allLabels.mat'], 'allAreaLabels', 'animalLabels', 'groupLabels')
sortT = NaN(1,length(allAreaLabels));

clustFile = dir([umapPath '20200619T2320\umap_area*.mat']);
clustFile = {clustFile(:).name};
checkFile = dir([umapPath '*checked.mat']);
checkFile = {checkFile(:).name};
h = figure('renderer','painters');
for iArea = unique(allAreaLabels)
    areaIdx = find(allAreaLabels == iArea);
    
%     cFile = [umapPath checkFile{contains(checkFile, ['area' int2str(iArea) '_checked.mat'])}]; %get file
    cFile = [umapPath '20200619T2320\umap_area_' num2str(iArea) '.mat'];
    load(cFile); disp(cFile);
    newT = removeGaps(spatialLabels);
    
    if isnan(max(sortT))
        sortT(areaIdx(cIdx)) = newT;
    else
        sortT(areaIdx(cIdx)) = max(sortT) + newT;
    end
    
    % show clusters
    if plotUMAP
        figure(h);
        cFile = [umapPath clustFile{contains(clustFile, ['area' int2str(iArea) '_clustered.mat'])}]; %get file
        load(cFile); disp(cFile);
        scatter(umapOut(hdbLabels==-1,1),umapOut(hdbLabels==-1,2),100,[0.75, 0.75, 0.75],'.'); axis square; hold on;
        scatter(umapOut(hdbLabels>-1,1),umapOut(hdbLabels>-1,2),100,hdbLabels(hdbLabels>-1),'.'); axis square
        title(iArea); colormap(jet(20)); hold off;
        pause;
    end
end
close(h);

%%
Cnt = 0;
clustIdx = []; clustCnt = []; clust = [];
allClust = unique(sortT); %all clusters IDs
allClust = allClust(~isnan(allClust));
for x = allClust
    Cnt = Cnt + 1;
    for iAnimals = unique(animalLabels)
        clust(iAnimals,Cnt) = sum(sortT(animalLabels == iAnimals) == x) ./ length(dimNr{iAnimals});
    end
    for iGroups = unique(groupLabels)
        cGroup = contains(animals, groups{iGroups});
        clustCnt(iGroups,Cnt) = mean(clust(cGroup,Cnt));
    end
    clustIdx(:,Cnt) = clustCnt(:,Cnt) ./ sum(clustCnt(:,Cnt));
end


%% sort groups based on clustering index
sortIdx = []; sortL = [];
for x = 1 : length(groups)
    [a,b] = sort(clustIdx(x,:),'descend');
    b = b(a > clustThresh);
    b(ismember(b,sortIdx)) = [];
    sortIdx = [sortIdx b];
    sortL(x) = length(sortIdx);
end
a = 1:size(clustIdx,2);
sortIdx = [sortIdx a(~ismember(a,sortIdx))];

figure('renderer','painters');
subplot(2,3,1:2); colormap viridis;
imagesc(clust(:,sortIdx)); 
% axis image;
nvline(sortL+0.5, 'w','linewidth',2);
nhline([4 6 8] + 0.5, 'w','linewidth',2);
caxis([0 1]); xlim([1 sortL(end)])

subplot(2,3,3);
imagesc(clust(:,sortIdx)); 
% axis image;
caxis([0 1]); xlim([sortL(end)+1 size(clust,2)])

subplot(2,3,4:5);
cData = clust ./ nanstd(clust, [], 1);
imagesc(cData(:,sortIdx)); 
% axis image;
nvline(sortL+0.5, 'w','linewidth',2);
nhline([4 6 8] + 0.5, 'w','linewidth',2);
caxis([0 3]); 
xlim([1 sortL(end)])

subplot(2,3,6);
imagesc(cData(:,sortIdx)); 
caxis([0 3]);
% axis image;
xlim([sortL(end)+1 size(clust,2)])

%% compare selective clusters and get averages
clustIdxAll = cell(1,length(groups));
clustDiff = cell(1,length(groups));
clustCorr = cell(1,length(groups));
clustMaps = cell(1,length(groups));
clustID = cell(1,length(groups));
for x = 1 : length(groups)
   
    cClust = find(clustIdx(x,:) > clustThresh); %most specific clusters for this group
    oClust = find(clustIdx(x,:) < clustThresh); %other clusters
    clustIdxAll{x} = clustIdx(x,cClust);
    clustID{x}(:,1) = allClust(cClust);
    
    % rebuild data
    cData = arrayfun(@(x) nanmean(X(:,sortT == allClust(x),:),2)', cClust, 'UniformOutput', false);
    cData = maxnorm(cat(1,cData{:}),2);
    oData = arrayfun(@(x) nanmean(X(:,sortT == allClust(x),:),2)', oClust, 'UniformOutput', false);
    oData = maxnorm(cat(1,oData{:}),2);
    
    % compare current regressors agains other regs
    for y = 1 : length(cClust)
        
%         % alternative metric to correlation. doesnt seem to add anything.
%         temp = cData(y,:) - oData;
%         temp = arrayShrink(temp',nShrinkMask,'split'); %recreate difference images to do smoothing
%         for xx = 1:size(temp,3)
%             temp(:,:,xx) = smoothImg(temp(:,:,xx), 'gaussian', 3, 1.76);
%         end
%         temp = arrayShrink(temp,nShrinkMask,'merge');
%         clustCorr{x}(y,:) = max(temp, [], 1) - min(temp, [], 1);
        
        temp = corrcoef(cat(1,cData(y,:),oData)');
        clustCorr{x}(y,:) = temp(1,2:end);
        clustMaps{x}{y,1} = cData(y,:)';
        
        % get 2 closest clusters
        [~,temp] = sort(clustCorr{x}(y,:),'descend');
        clustID{x}(y,2:3) = allClust(oClust(temp(1:2)));
        clustMaps{x}{y,2} = oData(temp(1:2),:)';
    end
    
    fprintf('done: group %d - %d / %d\n',x,y,length(cClust));
end


%% show comparison vs selectivity
figure('renderer','painters');
for x = 1 : 4
    subplot(1,4,x)
    plot(max(clustCorr{x},[],2),clustIdxAll{x},'o'); title(groups{x}); 
    xlabel('Similarity'); ylabel('Specificity'); axis square; niceFigure;
    xlim([0.4 1.1]); ylim([0.4 1.1]);hline(0.5)
end

%% show potential clusters - EMX
figure('renderer','painters');
iAnimal = 1;
temp = find(max(clustCorr{iAnimal},[],2)' <0.9 & clustIdxAll{iAnimal} > 0.8);
[~,b] = sort(max(clustCorr{iAnimal}(temp,:),[],2)', 'ascend');
temp = temp(b(1:2));
% temp(2) = []; %dont show split cluster
cCnt = length(temp);
bkgrnd = NaN(size(shrinkMask));
bkgrnd(~shrinkMask) = 0;
for x = 1 : cCnt
    
    subplot(cCnt,3,((x-1)*3)+1);
    bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),1};
    cImg = imageScale(bkgrnd,0.9); axis image
    cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
    title(['targetCluster: ' num2str(clustID{iAnimal}(temp(x))) ', Specificity: ' num2str(round(clustIdxAll{iAnimal}(temp(x)),2))]);
    
    for y = 1 : 2
        subplot(cCnt,3,((x-1)*3)+y+1);
        bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),2}(:,y);
        cImg = imageScale(bkgrnd,0.9); axis image
        cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
        title(round(corr2(clustMaps{iAnimal}{temp(x),1},clustMaps{iAnimal}{temp(x),2}(:,y)),2));
%         title(['nearCluster: ' num2str(clustID{iAnimal}(temp(x)))]);
    end
end
disp('Occurence probability for selected clusters:')
disp(clust(:,ismember(clustIdx(1,:), clustIdxAll{iAnimal}(temp))));

%% show potential clusters - Fez
figure('renderer','painters');
iAnimal = 2;
temp = find(max(clustCorr{iAnimal},[],2)' <0.9 & clustIdxAll{iAnimal} > 0.8);
[~,b] = sort(max(clustCorr{iAnimal}(temp,:),[],2)', 'ascend');
temp = temp(b(1:2));
% temp(2) = []; %dont show split cluster
cCnt = length(temp);
for x = 1 : cCnt
    
    subplot(cCnt,3,((x-1)*3)+1);
    bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),1};
    cImg = imageScale(bkgrnd,0.9); axis image
    cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
    title(['targetCluster: ' num2str(clustID{iAnimal}(temp(x))) ', Specificity: ' num2str(round(clustIdxAll{iAnimal}(temp(x)),2))]);
    
    for y = 1 : 2
        subplot(cCnt,3,((x-1)*3)+y+1);
        bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),2}(:,y);
        cImg = imageScale(bkgrnd,0.9); axis image
        cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
        title(round(corr2(clustMaps{iAnimal}{temp(x),1},clustMaps{iAnimal}{temp(x),2}(:,y)),2));
%         title(['nearCluster: ' num2str(clustID{iAnimal}(temp(x)))]);
    end
end
disp('Occurence probability for selected clusters:')
disp(clust(:,ismember(clustIdx(1,:), clustIdxAll{iAnimal}(temp))));

%% show potential clusters - Plexin
figure('renderer','painters');
iAnimal = 3;
temp = find(max(clustCorr{iAnimal},[],2)' <1 & clustIdxAll{iAnimal} > 0.8);
[~,b] = sort(max(clustCorr{iAnimal}(temp,:),[],2)', 'ascend');
temp = temp(b([2 5]));
% temp(2) = []; %dont show split cluster
cCnt = length(temp);
for x = 1 : cCnt
    
    subplot(cCnt,3,((x-1)*3)+1);
    bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),1};
    cImg = imageScale(bkgrnd,0.9); axis image
    cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
    title(['targetCluster: ' num2str(clustID{iAnimal}(temp(x))) ', Specificity: ' num2str(round(clustIdxAll{iAnimal}(temp(x)),2))]);
    
    for y = 1 : 2
        subplot(cCnt,3,((x-1)*3)+y+1);
        bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),2}(:,y);
        cImg = imageScale(bkgrnd,0.9); axis image
        cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
        title(round(corr2(clustMaps{iAnimal}{temp(x),1},clustMaps{iAnimal}{temp(x),2}(:,y)),2));
%         title(['nearCluster: ' num2str(clustID{iAnimal}(temp(x)))]);
    end
end
disp('Occurence probability for selected clusters:')
disp(clust(:,ismember(clustIdx(1,:), clustIdxAll{iAnimal}(temp))));

%% show potential clusters - CSP
figure('renderer','painters');
iAnimal = 4;
temp = find(max(clustCorr{iAnimal},[],2)' <0.9 & clustIdxAll{iAnimal} > 0.7);
[~,b] = sort(max(clustCorr{iAnimal}(temp,:),[],2)', 'ascend');
temp = temp(b([1 2]));
cCnt = length(temp);
bkgrnd = NaN(size(shrinkMask));
for x = 1 : cCnt
    
    subplot(cCnt,3,((x-1)*3)+1);
    bkgrnd = NaN(size(shrinkMask));
    bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),1};
    cImg = imageScale(bkgrnd,0.9); axis image
    cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
    title(['targetCluster: ' num2str(clustID{iAnimal}(temp(x))) ', Specificity: ' num2str(round(clustIdxAll{iAnimal}(temp(x)),2))]);
    
    for y = 1 : 2
        subplot(cCnt,3,((x-1)*3)+y+1);
        bkgrnd = NaN(size(shrinkMask));
        bkgrnd(~nShrinkMask) = clustMaps{iAnimal}{temp(x),2}(:,y);
        cImg = imageScale(bkgrnd,0.9); axis image
        cImg.Parent.CLim(1) = 0; colormap(inferno(256)); rateDisc_plotAllenOutline
        title(round(corr2(clustMaps{iAnimal}{temp(x),1},clustMaps{iAnimal}{temp(x),2}(:,y)),2));
%         title(['nearCluster: ' num2str(clustID{iAnimal}(temp(x)))]);
    end
end
