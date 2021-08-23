% function rateDisc_dimensionalityCompare
% code to assess the dimensionality of the widefield data for all animals
% of a given group. This uses rateDisc_plotAreaSVD, which
% looks at PCs when running SVD over all recordings across training. The
% output is a boxplot comparing dimensionality of EMX, FEZ and Plexin mice.

[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
lastRec = 10; %limit session to end of auditory discrimination
animals = dataOverview(1:10,1);
groups = {'mSM', 'Fez', 'Plex' 'CSP'};
groupColor = {[1 0 0] [0 0 1] [0 0 0] [0 1 0]}; %colors for groups
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
[xRange, yRange] = rateDisc_maskRange(allenMask); % get inner range of allenMask
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
cPath = 'Q:\BpodImager\Animals\'; %data path on the server
trainingRange = 'allAudio';

%% get data from all mice
Cnt = 0;
groupCnt = cell(1,length(groups));
gDimCnt = []; lDimCnt = []; recDates = []; recType = []; totalVar = []; aDetection = []; groupId = [];

for iAnimals = 1 : length(animals)
    cAnimal = animals{iAnimals}; % current animal

    % check if animal belongs to a group
    for iGroups = 1 : length(groups)
        if contains(cAnimal, groups{iGroups})
            [cgDimCnt, clDimCnt, sessionVar, sessionVarMap, sessionAvgVarMap, sessionModelVarMaps, recs] = rateDisc_plotAreaSVD(cPath, cAnimal, lastRec, false,trainingRange);
            [Performance,bhv] = rateDisc_sessionPerformance(cAnimal,cPath,recs);
            
            Cnt = Cnt + 1;
            if ~isempty(cgDimCnt)
                
                gDimCnt = [gDimCnt cgDimCnt];
                lDimCnt = [lDimCnt clDimCnt];
                totalVar = [totalVar sessionVar(end,:)];
                groupId = [groupId ones(1,length(cgDimCnt))*Cnt];
                groupCnt{iGroups} = [groupCnt{iGroups} Cnt];
                varMaps{Cnt} = sessionVarMap;
                varAvgMaps{Cnt} = sessionAvgVarMap;
                modelVarMaps{Cnt} = sessionModelVarMaps;
                
                recDates = [recDates Performance.date];
                aDetection = [aDetection Performance.Detection(1,:)];
                recType = [recType Performance.recType];
                
            end
        end
    end
end

%% variance figure
gPos = zeros(1,length(varMaps));
for iGroups = 1 : length(groupCnt)
    gPos(groupCnt{iGroups}) = 1 + max(gPos) + (1:length(groupCnt{iGroups}));
end
gPos(gPos == 0) = [];
    
figure('renderer','painters');
subplot(1,3,1);
boxplot(totalVar, groupId, 'positions', gPos, 'labels', animals(unique(groupId)))
axis square; ylim([0 75]); title('variance compare');
ylabel('total variance'); niceFigure(gca);

subplot(1,3,2); hold on
xMax = 0;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
        cDate = recDates(groupId == groupCnt{iGroups}(iAnimals));
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
        plot(cDate, smooth(totalVar(groupId == groupCnt{iGroups}(iAnimals)),3), 'Color', cColor);
    end
end
axis square; ylabel('session variance'); title('variance/sessions');
xlabel('training days'); niceFigure(gca); ylim([0 60]); xlim([0 xMax]);

subplot(1,3,3); hold on
xMax = 0;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
        cDate = recDates(groupId == groupCnt{iGroups}(iAnimals));
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
        plot(cDate, smooth(aDetection(groupId == groupCnt{iGroups}(iAnimals)),3), 'Color', cColor);
    end
end
axis square; ylabel('detection performance'); title('performance/sessions');
xlabel('training days'); niceFigure(gca); ylim([0.4 1]); xlim([0 xMax]);

%% learning figure
figure('renderer','painters');
xMax = 0;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
        subplot(1,2,1); hold on
        cIdx = groupId == groupCnt{iGroups}(iAnimals) & recType < 7;
        cIdx = cIdx & (recDates - min(recDates(cIdx))) < 30;
        
        cDate = recDates(cIdx);
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
        plot(cDate, smooth(aDetection(cIdx),3), 'Color', cColor);
        axis square; ylabel('detection performance'); title('detection/learning');
        xlabel('training days'); niceFigure(gca); ylim([0.4 1]); xlim([0 xMax]);
        
        subplot(1,2,2); hold on
        cIdx = groupId == groupCnt{iGroups}(iAnimals) & recType == 5;
        cIdx = cIdx & (recDates - min(recDates(cIdx))) <= 20;
        cDate = recDates(cIdx);
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
        plot(cDate, smooth(aDetection(cIdx),3), 'Color', cColor);
        axis square; ylabel('detection performance'); title('discrimination');
        xlabel('training days'); niceFigure(gca); ylim([0.4 1]); xlim([0 xMax]);
    end
end

%% show variance difference
flipHS = false; %flag to collapse both hemisspehres in one
segSplit = 3; %split sequences after this recType (3 is 90th prctile)
segL = [10 10]; %number of sessions before and after split
colormap(inferno(256));
segT = {'Before' 'After'};
segL = [10 2 10];
circleR = 30; %radius of circle mask
circleP = [215 310]; %position of circle mask
% circleP = [80 380]; %position of circle mask
% circleP = [240 150]; %position of circle mask
if flipHS
    cMask = allenMask(:,1:round(size(allenMask,2)/2));
else
    cMask = allenMask;
end
cMask = ~createCirclesMask(size(cMask),circleP,circleR);

for iGroups = 1:4
    figure('renderer','painters'); hold on
    cData = nan([size(cMask) 2]);
    for iSegs = 1 : 2
        for iAnimals = 1 : length(groupCnt{iGroups})
            if iSegs == 1
                cRange = find(recType(groupId == groupCnt{iGroups}(iAnimals)) <= segSplit);
            else
                cRange = find(recType(groupId == groupCnt{iGroups}(iAnimals)) > segSplit);
            end
            cRange = cRange(end-min([segL(1) length(cRange)])+1:end);
            temp = rateDisc_removeOutline(varMaps{groupCnt{iGroups}(iAnimals)}(:,:,cRange),10,0);
            temp = nanmean(temp,3);
            if flipHS
                temp = rateDisc_flipHS(temp,'mean');
            end
            cData(:,:,iSegs) = nanmean(cat(3,cData(:,:,iSegs),temp),3);
        end
        subplot(1,4,iSegs);
        cImg = imageScale(cData(:,:,iSegs)); axis image;
        temp = cData(:,:,1); caxis([0 prctile(temp(:),95)]);
        colormap(cImg.Parent,inferno(256))
        title([groups{iGroups} '-' segT{iSegs}])
    end
    subplot(1,4,iSegs+1);
    cImg = imageScale(diff(cData,[],3)./cData(:,:,1),0.6); axis image;
    colormap(cImg.Parent,colormap_blueblackred(256))
    title([groups{iGroups} '- Percent change'])
    hold on;
    rectangle('Position',[circleP-circleR [circleR circleR].*2],'Curvature',1,'EdgeColor','w');
    niceFigure;

    % get activity in ROI
    cData = NaN(length(groupCnt{iGroups}),sum(segL)+1);
    for x = 1 : length(groupCnt{iGroups})
        if flipHS
            cVar = rateDisc_flipHS(varMaps{groupCnt{iGroups}(x)},'mean');
        else
            cVar = varMaps{groupCnt{iGroups}(x)};
        end
        cVar = nanmean(arrayShrink(cVar,cMask,'merge'),1); %variance in area of interest
        cRecs = recType(groupId == groupCnt{iGroups}(x)); %rectypes for current animal

        cRange = find(cRecs <= segSplit);
        cRange = cRange(end-min([segL(1) length(cRange)])+1:end); %only use latest recs
        cData(x, segL(1)-length(cRange)+1:segL(1)) = cVar(cRange);
        
        cRange = find(cRecs <= segSplit + 1 & cRecs > segSplit);
        cRange = cRange(end-min([segL(2) length(cRange)])+1:end);
        cData(x, segL(1)+1:segL(1) + length(cRange)) = cVar(cRange);
        
        cRange = find(cRecs > segSplit+1);
        cRange = cRange(end-min([segL(3) length(cRange)])+1:end);
        cData(x, sum(segL(1:2))+2:sum(segL(1:2)) + 1 + length(cRange)) = cVar(cRange);
    end
    subplot(1,4,4);
    cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
    plot(1:length(cData),cData, 'Color', cColor); axis square; hold on;
    plot(1:length(cData),mean(cData,1), 'Color', cColor, 'lineWidth',4);
    title('Variance change - detection/discrimination');
    ylim([0 1.2E-3])
end
    
%% show modeled variance maps
h = figure('renderer', 'painters');
cReg = 5; %discrimination only
for iGroups = 1 : length(groupCnt)
    for x = 1 : 3
        subplot(length(groupCnt), 3, (iGroups-1)*3 + x);
        
        Cnt = 0;
        cMap = NaN(size(allenMask,1), size(allenMask,2), length(groupCnt{iGroups}));
        avgVar = cell(1, length(groupCnt{iGroups}));
        for iAnimals = groupCnt{iGroups}
            Cnt = Cnt + 1;            
            avgVar{Cnt} = squeeze(nanmean(arrayShrink(modelVarMaps{iAnimals}(:,:,x,:), allenMask, 'merge'),1));
            cIdx = recType(groupId == groupCnt{iGroups}(Cnt));
            cIdx = cIdx == cReg;
            if x == 1
                cMap(:,:,Cnt) = nanmean(modelVarMaps{iAnimals}(:,:,1,cIdx),4);
            else
                cMap(:,:,Cnt) = nanmean(modelVarMaps{iAnimals}(:,:,1,cIdx),4) - nanmean(modelVarMaps{iAnimals}(:,:,x,cIdx),4);
            end
        end
        allVar{x, iGroups} = cat(1,avgVar{:});
        cMap = nanmean(cMap,3);
        
        cRange = [0 prctile(cMap(:),99)];
        mapImg = imshow(cMap,cRange);
        colormap(inferno); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(groups{iGroups}, 'FontSize', 14)
        
        rateDisc_plotAllenOutline(gca,'L');
    end
end

%% model variance boxplot - full, movement and task
clear cData
newId = [];
Cnt = 1;
for x = 1 : length(groups)
    newId = [newId ones(1,size(allVar{1,x},1)).*Cnt];
    Cnt = max(newId(:)) + 1; disp(Cnt);
end
 
% cData(:,1) = cat(1,allVar{1,:});
% cData(:,2) = cat(1,allVar{3,:}); %move
% cData(:,3) = cat(1,allVar{2,:}); %task
% cData(:,2) = cat(1,allVar{1,:}) - cat(1,allVar{2,:}); %unique move
% cData(:,3) = cat(1,allVar{1,:}) - cat(1,allVar{3,:}); %unique task

h = figure('renderer', 'painters'); hold on
cLabels = {'Full model','Movement', 'Task'};
for x = 1 : 4
    
    cIdx = recType(newId == x) > 2;
    cData = cat(2, allVar{:,x});
    cData = reshape(cData(cIdx,[1 3 2]),[],1);
    cID = reshape(ones(1,sum(cIdx))' .* (1:3), [], 1);
    
    subplot(1,4,x);
    violinplot(cData,cID); axis square
    ax = gca; ax.XTickLabel = cLabels; ax.TickLength = [0 0];
    ylabel('cvR^2'); title(groups{x}); ax.YLim = [0 0.8];
    grid on; 
end

% %same as above but for different models grouped together
% h = figure('renderer', 'painters'); hold on
% cLabels = {'Full model','Movement', 'Task'};
% for x = 1 : 3
%     
%     cData = cat(1, allVar{x,:});
%     cIdx = recType > 2;
% 
%     subplot(1,3,x);
%     violinplot(cData(cIdx),newId(cIdx)); axis square
%     ax = gca; ax.XTickLabel = cLabels; ax.TickLength = [0 0];
%     ylabel('cvR^2'); title(groups{x}); ax.YLim = [0 0.8];
%     grid on; 
% end

% boxplot(cat(1,allVar{1,:}) - cat(1,allVar{2,:}), groupId, 'Color', 'b', 'positions', gPos, 'labels', animals(unique(groupId)))  %unique motor 
% boxplot(cat(1,allVar{1,:}) - cat(1,allVar{3,:}), groupId, 'Color', 'k', 'positions', gPos, 'labels', animals(unique(groupId)))  %unique task

%% split maps by training progress
segSplit = 4; %split sequences after this recType (3 is 90th prctile)
segL = [5 5 5];

h = figure('renderer', 'painters');
for iGroups = 1:4
    for iSegs = 2:3
        cMap = NaN(size(allenMask,1),size(allenMask,2),length(groupCnt{iGroups}));
        for iAnimals = 1 : length(groupCnt{iGroups})
            cRecs = recType(groupId == groupCnt{iGroups}(iAnimals));
            if iSegs == 1
%                 cRange = find(cRecs <= segSplit);
%                 cRange = cRange(1:min([segL(iSegs) length(cRange)]));
                cRange = 1:length(segL(1));
%                 cRange = find(cRecs <= segSplit);
%                 cRange = cRange(1:min([segL(iSegs) length(cRange)]));
%                 cRange = cRange(end-min([segL(iSegs) length(cRange)])+1:end);
            elseif iSegs == 2
%                 cRange = find(cRecs >= segSplit & cRecs <= segSplit+2);
%                 cRange = cRange(end-min([segL(iSegs) length(cRange)])+1:end);
%                 cRange = cRange(1:min([segL(iSegs) length(cRange)]));
                cRange = find(cRecs <= segSplit);
%                 cRange = find(ismember(cRecs, [3]));
                cRange = cRange(end-min([segL(iSegs) length(cRange)])+1:end);
            elseif iSegs == 3
%                 cRange = find(cRecs > segSplit+1);
%                 cRange = cRange(end-min([segL(iSegs) length(cRange)])+1:end);
%                 cRange = cRange(1:min([segL(iSegs) length(cRange)]));
                cRange = find(ismember(cRecs, segSplit+1));
                cRange = cRange(end-min([segL(iSegs) length(cRange)])+1:end);

            end
%             cMap(:,:,iAnimals) = nanmean(modelVarMaps{groupCnt{iGroups}(iAnimals)}(:,:,2,cRange),4);
            cMap(:,:,iAnimals) = nanmean(modelVarMaps{groupCnt{iGroups}(iAnimals)}(:,:,1,cRange) - modelVarMaps{groupCnt{iGroups}(iAnimals)}(:,:,3,cRange),4);
        end
        
        %
        subplot(2, length(groupCnt), (iGroups-1)*2 + iSegs - 1);
        cImg = imageScale(nanmean(cMap,3)); cImg.Parent.CLim(1)=0;
        rateDisc_plotAllenOutline(gca); title(groups{iGroups});
%         caxis([0 0.5]);
    end
end
colormap(inferno(256));
% colormap(colormap_blueblackred(256));


%%
figure('renderer','painters'); hold on


%% show variance in target area (split detection and discrimination)
figure('renderer','painters'); hold on
segSplit = 3; %split sequences after this recType (3 is 90th prctile)
segL = [10 2 10];
for iSegs = 1:3
    for iGroups = 1:4
        for iAnimals = 1 : length(groupCnt{iGroups})
            cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
            
            cData = NaN(1,segL(iSegs));
            if iSegs == 1
                cRange = find(recType <= segSplit & groupId == groupCnt{iGroups}(iAnimals));
                cData(end-min([segL(iSegs) length(cRange)])+1:end) = totalVar(cRange(end-min([segL(iSegs) length(cRange)])+1:end));
                cRange = find(recType > segSplit & recType <= segSplit+1 & groupId == groupCnt{iGroups}(iAnimals));
                cData = [cData totalVar(cRange(1:min([segL(iSegs+1) length(cRange)])))];
                plot(1:length(cData),cData, 'Color', cColor);
%             elseif iSegs == 2
%                 cRange = find(recType > segSplit & recType <= segSplit+1 & groupId == groupCnt{iGroups}(iAnimals));
%                 cData(1:min([segL(iSegs) length(cRange)])) = totalVar(cRange(1:min([segL(iSegs) length(cRange)])));
%                 plot(segL(1)+1:sum(segL(1:2)),cData, 'Color', cColor);
            elseif iSegs == 3
                cRange = find(recType > segSplit + 1 & groupId == groupCnt{iGroups}(iAnimals));
                cData(1:min([segL(iSegs) length(cRange)])) = totalVar(cRange(1:min([segL(iSegs) length(cRange)])));
                plot(sum(segL(1:2))+1:sum(segL(1:3)),cData, 'Color', cColor);
            end
        end
        
    end
    axis square; ylabel('session variance'); title('variance/sessions');
    niceFigure(gca); 
%     ylim([0 75]);
end


%% show average variance map for each animal - averaged over all sessions
figure('renderer', 'painters');
% bIdx{1} = [4 5 6 7 8 9 10 12 14 15 16 19]; %csp22
% bIdx{2} = [5 9 10 12 13 14 16]; %csp23
Cnt = 0;
for iAnimals = [groupCnt{:}]
    Cnt = Cnt + 1;
    subplot(2, ceil(length([groupCnt{:}])/2), Cnt);
    
    cIdx = 1 : size(varMaps{Cnt},3);
%     cIdx(bIdx{iAnimals}) = [];
    cMap = nanmean(varMaps{Cnt}(:,:,cIdx),3);
    cRange = [0 prctile(cMap(:),99)];
    mapImg = imshow(cMap,cRange);
    
    colormap(inferno); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(animals{iAnimals}, 'FontSize', 14)
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% group average variance
figure('renderer', 'painters');
for iGroups = 1 : length(groupCnt)
    subplot(2, 2, iGroups);

    Cnt = 0;
    cMap = NaN(size(allenMask,1), size(allenMask,2), length(groupCnt{iGroups}));
    for iAnimals = groupCnt{iGroups}
        Cnt = Cnt + 1;
        cMap(:,:,Cnt) = nanmean(varMaps{iAnimals},3);
    end
    cMap = nanmean(cMap,3);
    
    cRange = [0 prctile(cMap(:),99)];
    mapImg = imshow(cMap,cRange);
    colormap(inferno); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(animals{iAnimals}, 'FontSize', 14)
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

% subplot(2, 2, iGroups+1); hold on
% for iGroups = 1 : length(groupCnt)
%     Cnt = 0;
%     cData = NaN(200, length(groupCnt{iGroups}));
%     for iAnimals = groupCnt{iGroups}
%         Cnt = Cnt + 1;
%         cDate = recDates(groupId == iAnimals);
%         cDate = cDate - min(cDate) + 1;
%         cData = totalVar(groupId == groupCnt{iGroups}(Cnt));
%         
%         cData(1:cDate(end),Cnt) = interp1(cDate,cData,1:cDate(end),'linear');
%     end
%     cLines(iGroups) = stdshade(cData(~isnan(nanmean(cData,2)), :)', find(~isnan(nanmean(cData,2))), groupColor{iGroups}, 0.5, 3);
% end
% groups{1} = 'Emx';
% title('Session variance per group');
% axis square; legend(cLines, groups);
% ylabel('session variance'); 
% xlabel('training days'); niceFigure(gca);

%% group average variance change
figure('renderer', 'painters');
for iGroups = 1 : length(groupCnt)
    subplot(1, length(groupCnt), iGroups);

    Cnt = 0;
    cMap = NaN(size(allenMask,1), size(allenMask,2), length(groupCnt{iGroups}));
    for iAnimals = groupCnt{iGroups}
        Cnt = Cnt + 1;
        cMap(:,:,Cnt) = nanvar(varMaps{iAnimals},[],3);
    end
    cMap = nanmean(cMap,3);
    
    cRange = [0 prctile(cMap(:),95)];
    mapImg = imshow(cMap,cRange);
    colormap(inferno); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(animals{iAnimals}, 'FontSize', 14)
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% dimensionality figure
gPos = zeros(1,length(groupCnt));
for iGroups = 1 : length(groupCnt)
    gPos(groupCnt{iGroups}) = 1 + max(gPos) + (1:length(groupCnt{iGroups}));
end
gPos(gPos == 0) = [];
    
figure;
subplot(2,3,1); hold on;
boxplot(gDimCnt, groupId, 'positions', gPos, 'labels', animals(unique(groupId)))
ylim([0 250]); title('across session-PCs');
ylabel('nr of PCs'); niceFigure(gca); axis square; 


subplot(2,3,2); hold on;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
        cDate = recDates(groupId == groupCnt{iGroups}(iAnimals));
        plot(cDate - min(cDate) + 1, smooth(gDimCnt(groupId == groupCnt{iGroups}(iAnimals)),3), 'Color', groupColor{iGroups});
    end
end
axis square; ylabel('nr of PCs'); title('global PCs for each session');
xlabel('training days'); niceFigure(gca);


subplot(2,3,2); hold on;
xMax = 0;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
        cDate = recDates(groupId == groupCnt{iGroups}(iAnimals));
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        plot(cDate, smooth(gDimCnt(groupId == groupCnt{iGroups}(iAnimals)),3), 'Color', groupColor{iGroups});
    end
end
axis square; ylabel('nr of PCs'); xlim([0 xMax]);
xlabel('training days'); niceFigure(gca);  title('global PCs for each session');

subplot(2,3,3); hold on;
for iGroups = 1 : length(groupCnt)
    Cnt = 0;
    cData = NaN(200, length(groupCnt{iGroups}));
    for iAnimals = groupCnt{iGroups}
        Cnt = Cnt + 1;
        cDate = recDates(groupId == iAnimals);
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        cData(1:cDate(end),Cnt) = interp1(cDate,gDimCnt(groupId == groupCnt{iGroups}(Cnt)),1:cDate(end),'linear');
    end
    cLines(iGroups) = stdshade(cData(~isnan(nanmean(cData,2)), :)', find(~isnan(nanmean(cData,2))), groupColor{iGroups}, 0.5, 15);
end
axis square; ylabel('nr of PCs'); xlim([0 xMax]);
xlabel('training days'); niceFigure(gca);

subplot(2,3,4); hold on;
boxplot(lDimCnt, groupId, 'positions', gPos, 'labels', animals(unique(groupId)))
ylim([0 250]);  title('within session-PCs');
ylabel('nr of PCs'); niceFigure(gca); axis square; 


subplot(2,3,5); hold on;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
        cDate = recDates(groupId == groupCnt{iGroups}(iAnimals));
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        plot(cDate - min(cDate) + 1, smooth(lDimCnt(groupId == groupCnt{iGroups}(iAnimals)),1), 'Color', groupColor{iGroups});
    end
end
axis square; ylabel('nr of PCs'); title('local PCs for each session'); xlim([0 xMax]);
xlabel('training days'); niceFigure(gca);


subplot(2,3,6); hold on;
for iGroups = 1 : length(groupCnt)
    Cnt = 0;
    cData = NaN(200, length(groupCnt{iGroups}));
    for iAnimals = groupCnt{iGroups}
        Cnt = Cnt + 1;
        cDate = recDates(groupId == iAnimals);
        cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
        cData(1:cDate(end),Cnt) = interp1(cDate,lDimCnt(groupId == groupCnt{iGroups}(Cnt)),1:cDate(end),'linear');
    end
    cLines(iGroups) = stdshade(cData(~isnan(nanmean(cData,2)), :)', find(~isnan(nanmean(cData,2))), groupColor{iGroups}, 0.5, 15);
end
axis square; ylabel('nr of PCs'); xlim([0 xMax]);
xlabel('training days'); niceFigure(gca); title('local PCs for each session');