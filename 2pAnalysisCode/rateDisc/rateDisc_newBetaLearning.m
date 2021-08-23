% function rateDisc_newBetaLearning
% code to assess the beta kernels for all animals of a given group. This
% one is for the paper and should show change in beta kernels and average
% over sessions. e.g. for choice and stimulus.

[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
lastRec = 5; %limit session to end of auditory discrimination
animals = dataOverview(:,1);
groups = {'mSM', 'Fez', 'Plex'};
groupColor = {[1 0 0] [1 0 1] [0 1 0]}; %colors for groups
groupCnt = cell(1,length(groups));
avgBeta = cell(1, length(groups));
load('allenDorsalMapSM.mat');
allenMask = dorsalMaps.allenMask;
load('trimMap.mat');
trimMap = arrayShrink(trimMap, allenMask, 'merge');
allenMask = dorsalMaps.allenMask;
[xRange, yRange] = rateDisc_maskRange(allenMask); % get inner range of allenMask
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
fileExt = 'org'; %non-orthogonalized model
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server

%%
for iAnimals = 1 : 8
    % for iAnimals = 1 : length(animals)
    cAnimal = animals{iAnimals}; % current animal
    bPath = [cPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, 'rateDisc'); %this sorts recordings by date

for iRecs = 1 : size(recs,1)
    % get some data
    load([bPath fileExt 'dimBeta.mat'])
    load([fPath 'regData.mat'], 'recIdx', 'idx', 'recLabels')
    load([fPath 'Vc.mat'], 'U')
    load([fPath 'opts2.mat'], 'opts')
    
    load([bPath fileExt 'betaMaps.mat'],'regBeta', 'regMeans', 'recs', 'regLabels');
    
    %remove sessions that are not present in 'varRecs'
    useIdx = false(1,length(recs));
    for iRecs = 1:length(recs)
       if any(ismember(cellstr(cat(1, varRecs.name)), recs(iRecs).name))
            useIdx(iRecs) = true;
       end
    end
    recs = recs(useIdx);
    regMeans = regMeans(:, :, useIdx);
    sessionBeta{iAnimals} = regMeans;

    fullMaps{iAnimals} = NaN(sum(~allenMask(:)), length(recs), 'single');
    motorMaps{iAnimals} = NaN(sum(~allenMask(:)), length(recs), 'single');
    taskMaps{iAnimals} = NaN(sum(~allenMask(:)), length(recs), 'single');
    for iGroups = 1 : length(groups)
        if contains(cAnimal, groups{iGroups})
            groupCnt{iGroups} = [groupCnt{iGroups} iAnimals];
            
            for iRegs = 1 : length(regBeta)
                if length(avgBeta{iGroups}) < iRegs
                    avgBeta{iGroups}{iRegs} = regBeta{iRegs};
                elseif size(regBeta{iRegs},2) > 0
                    if size(regBeta{iRegs},2) < size(avgBeta{iGroups}{iRegs},2)
                        regBeta{iRegs} = [regBeta{iRegs} zeros(sum(~allenMask(:)), size(avgBeta{iGroups}{iRegs},2) - size(regBeta{iRegs},2), 'single')]; %add columns for new regressors if neede
                    else
                        avgBeta{iGroups}{iRegs} = [avgBeta{iGroups}{iRegs} zeros(sum(~allenMask(:)),size(regBeta{iRegs},2) - size(avgBeta{iGroups}{iRegs},2), 'single')]; %add columns for new regressors if neede
                    end
                    avgBeta{iGroups}{iRegs} = avgBeta{iGroups}{iRegs} + ((regBeta{iRegs} - avgBeta{iGroups}{iRegs})/length(groupCnt{iGroups})); % update mean for current regressors
                end
            end
        end
    end
    
    % get performance for each recording
    Performance{iAnimals} = rateDisc_sessionPerformance(cAnimal,cPath,recs);
    
    % get model R2 for each recording
    for iRecs = 1 : length(recs)
        try
            fPath = [cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep];
            load([fPath 'mask.mat'], 'mask');
            load([fPath 'opts2.mat'], 'opts');
            
            load([fPath 'fullcorr.mat'], 'fullMap');
            fullMap = arrayShrink(fullMap, mask, 'split');
            fullMap = alignAllenTransIm(fullMap,opts.transParams); %align to allen
            fullMaps{iAnimals}(:,iRecs) = arrayShrink(fullMap.^2, allenMask, 'merge');
            
            load([fPath filesep 'motorregData.mat'], 'motorMap');
            motorMap = arrayShrink(motorMap, mask, 'split');
            motorMap = alignAllenTransIm(motorMap,opts.transParams); %align to allen
            motorMaps{iAnimals}(:,iRecs) = arrayShrink(motorMap.^2, allenMask, 'merge');
            
            load([fPath filesep 'taskregData.mat'], 'taskMap');
            taskMap = arrayShrink(taskMap, mask, 'split');
            taskMap = alignAllenTransIm(taskMap,opts.transParams); %align to allen
            taskMaps{iAnimals}(:,iRecs) = arrayShrink(taskMap.^2, allenMask, 'merge');
        end
    end
end


%% model reconstruction figure
figure('renderer','painters');
for iGroups = 1 : length(groups)
    
    fullMap = nanmean(cat(2, fullMaps{contains(animals, groups{iGroups})}),2);
    taskMap = nanmean(cat(2, taskMaps{contains(animals, groups{iGroups})}),2);
    motorMap = nanmean(cat(2, motorMaps{contains(animals, groups{iGroups})}),2);
    
    subplot(length(groups),3,iGroups);
    imageScale(arrayShrink(fullMap, allenMask, 'split'));
    colormap(inferno(256)); colorbar; caxis([0 0.5]);
    title([groups{iGroups} ' - full'], 'FontSize', 14)
    
    subplot(length(groups),3,iGroups + length(groups));
    imageScale(arrayShrink(fullMap - taskMap, allenMask, 'split'));
    colormap(inferno(256)); colorbar; caxis([0 0.4]);
    title([groups{iGroups} ' - unique motor'], 'FontSize', 14)
    
    subplot(length(groups),3,iGroups + (2 * length(groups)));
    imageScale(arrayShrink(fullMap - motorMap, allenMask, 'split'));
    colormap(inferno(256)); colorbar; caxis([0 0.01]);
    title([groups{iGroups} ' - unique task'], 'FontSize', 14)
    
end

%% Task model change figure
figure('renderer','painters');
for iGroups = 1 : length(groups)
    
    Cnt = 0;
    for iAnimals = groupCnt{iGroups}
        Cnt = Cnt + 1;
        cData = smooth(Performance{iAnimals}.Detection(1,:));
        minP = mean(cData(1:5));
        maxP = mean(cData(end-5:end));
        [a, b] = min(abs(cData - (((maxP - minP) / 2) + minP))); %find halfmax
        
        minData(:,Cnt) = nanmean(fullMaps{iAnimals}(:,max([1 b-7]):b-3),2) - nanmean(motorMaps{iAnimals}(:,max([1 b-7]):b-3),2);
        midData(:,Cnt) = nanmean(fullMaps{iAnimals}(:,b-2:b+2),2) - nanmean(motorMaps{iAnimals}(:,b-2:b+2),2);
        maxData(:,Cnt) = nanmean(fullMaps{iAnimals}(:,b+3:b+7),2) - nanmean(motorMaps{iAnimals}(:,b+3:b+7),2);
    end
    
    subplot(length(groups),3,iGroups);
    imageScale(arrayShrink(mean(minData,2), allenMask, 'split'));
    colormap(inferno(256)); colorbar; caxis([0 0.01]);
    title([groups{iGroups} ' - novice'], 'FontSize', 14)
    
    subplot(length(groups),3,iGroups + length(groups));
    imageScale(arrayShrink(mean(midData,2), allenMask, 'split'));
    colormap(inferno(256)); colorbar; caxis([0 0.01]);
    title([groups{iGroups} ' - learning'], 'FontSize', 14)
    
    subplot(length(groups),3,iGroups + (2 * length(groups)));
    imageScale(arrayShrink(mean(maxData,2), allenMask, 'split'));
    colormap(inferno(256)); colorbar; caxis([0 0.01]);
    title([groups{iGroups} ' - expert'], 'FontSize', 14)
    
end

%% motor/task maps across sessions
figure('renderer','painters');
Cnt = 0; xMax = 0;
for iGroups = 1 : length(groupCnt)
    for iAnimals = 1 : length(groupCnt{iGroups})
   
    
    cDate = Performance{groupCnt{iGroups}(iAnimals)}.date;
    cDate = cDate - min(cDate) + 1; xMax = max([xMax cDate(end)]);
    cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);
    
    subplot(1,2,1); hold on;
    cData = nanmean(fullMaps{groupCnt{iGroups}(iAnimals)},1) - nanmean(motorMaps{groupCnt{iGroups}(iAnimals)},1);
    cData = cData .* allVar{groupCnt{iGroups}(iAnimals)}(end,:);
    cData(cData > 1 | cData < -0.5) = nanmean(cData);
    plot(cDate, smooth(cData,3), 'Color', cColor);
    
    subplot(1,2,2); hold on;
    cData = nanmean(fullMaps{groupCnt{iGroups}(iAnimals)},1) - nanmean(taskMaps{groupCnt{iGroups}(iAnimals)},1);
    cData = cData .* allVar{groupCnt{iGroups}(iAnimals)}(end,:);
    plot(cDate, smooth(cData,3), 'Color', cColor);

    end
end
subplot(1,2,1); axis square; vline(0); title('Total task-explained variance'); 
xlabel('recordings'); ylabel('variance'); niceFigure(gca); 
subplot(1,2,2); axis square; vline(0); title('Total movement-explained variance'); 
xlabel('recordings'); ylabel('variance'); niceFigure(gca);


%% look at beta kernels
figure('renderer','painters');
cKernel = 'firstAudio';
% cKernel = 'firstAudio';
if strcmpi(cKernel, 'firstAudio')
    cKernel = {'lfirstAudStim' 'rfirstAudStim'};
elseif strcmpi(cKernel, 'licks')
    cKernel = {'lLick' 'rLick'};
elseif strcmpi(cKernel, 'handles')
    cKernel = {'lGrab' 'rGrab'};
end
rIdx = find(ismember(regLabels, cKernel));
cTrim = arrayShrink(trimMap, allenMask, 'split');

for iGroups = 1 : length(groupCnt)
    % find shortest kernel first if combining multiple ones
    minDur = inf; cMovie = [];
    for iRegs = rIdx
        minDur = min([minDur size(avgBeta{iGroups}{iRegs}, 2)]);
        if iRegs == rIdx(1)
            cMovie = avgBeta{iGroups}{iRegs};
        else
            cMovie = nansum(cat(3,cMovie(:, 1:minDur), avgBeta{iGroups}{iRegs}(:, 1:minDur)),3);
        end
    end
    %     compareMovie(arrayShrink(cMovie, allenMask,'split'), 'outline', dorsalMaps.edgeOutlineSplit)
    idx1 = 10:20;
    subplot(3,3, iGroups); hold on;
    [~, cRange] = imageScale(nanmean(arrayShrink(cMovie(:,idx1), allenMask,'split'),3), 99);
    colormap(colormap_blueblackred(256)); colorbar; title(groups{iGroups});
    
    idx2 = 40:50;
    subplot(3,3, iGroups + length(groups)); hold on;
    imageScale(nanmean(arrayShrink(cMovie(:,idx2), allenMask,'split'),3), 95);
    caxis(cRange); colormap(colormap_blueblackred(256)); colorbar;
%     
    % keep first response map for thresholding
    a = nanmean(arrayShrink(cMovie(:,idx1), allenMask,'split'),3);
    b = nanmean(arrayShrink(cMovie(:,idx2), allenMask,'split'),3);
    firstMap(:,:,iGroups) = a(:, 1:size(a,2)/2);
    secondMap(:,:,iGroups) = b(:, 1:size(b,2)/2);
end

% threshold maps to get some areas
cFrame = imadjust(firstMap(:,:,3)); %improve contrast
cFrame = cFrame > 0.80;
cFrame = bwareaopen(cFrame,100); %don't use areas that are smaller as 10 pixels
cFrame = imdilate(cFrame,strel('disk',6));
cFrame = imerode(cFrame,strel('disk',6));
cFrame(allenMask(:, 1:size(a,2)/2)) = false;
cFrame = bwareaopen(cFrame,300); %don't use areas that are smaller as 10 pixels
cAreas = outlineAndSmooth(cFrame);

cFrame = imadjust(firstMap(:,:,2)); %improve contrast
cFrame = cFrame > 0.85;
cFrame = bwareaopen(cFrame,100); %don't use areas that are smaller as 10 pixels
cFrame = imdilate(cFrame,strel('disk',6));
cFrame = imerode(cFrame,strel('disk',6));
cFrame(allenMask(:, 1:size(a,2)/2)) = false;
cFrame = bwareaopen(cFrame,400); %don't use areas that are smaller as 10 pixels
cAreas = [cAreas{1};outlineAndSmooth(cFrame)];

% plot areas outlines
for iAreas = 1 : length(cAreas)
   for x = 1 : 6
       subplot(3,3,x);
       plot(cAreas{iAreas}(:,2),cAreas{iAreas}(:,1), 'linewidth', 2);
   end
end

for iGroups = 1 : length(groupCnt)
    % find shortest kernel first if combining multiple ones
    minDur = inf; cMovie = [];
    for iRegs = rIdx
        minDur = min([minDur size(avgBeta{iGroups}{iRegs}, 2)]);
        if iRegs == rIdx(1)
            cMovie = avgBeta{iGroups}{iRegs};
        else
            cMovie = nansum(cat(3,cMovie(:, 1:minDur), avgBeta{iGroups}{iRegs}(:, 1:minDur)),3);
        end
    end
    cMovie = arrayShrink(cMovie, allenMask,'split');

    subplot(3,3, iGroups + length(groups)*2); hold on;
    for iAreas = 1 : length(cAreas)
        cMask = ~poly2mask(cAreas{iAreas}(:,2),cAreas{iAreas}(:,1),size(cFrame,1),size(cFrame,2));
        cTrace = nanmean(arrayShrink(cMovie(:, 1:size(cMovie,2)/2, :), cMask, 'merge'),1);
        plot(cTrace);
    end
    niceFigure(gca); axis square;
    vline([idx1(1) idx1(end) idx2(1) idx2(end)]);
end
