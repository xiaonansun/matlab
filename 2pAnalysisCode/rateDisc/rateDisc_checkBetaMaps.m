if ~exist('fileExt','var') || isempty(fileExt)
    fileExt = 'org'; %non-orthogonalized model by default
end

% cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
% cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
cPath = 'Q:\BpodImager\Animals\'; %path to churchlandNAS
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
varThresh = 4; %threshold in SDUs for recordings with too high variance (indicative of bad hemo correction)
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
trainingRange = 'audioDisc'; %use this to run analysis only in a certain range of training
leftHS = find(ismember(dorsalMaps.sidesSplit,'L')); %index for labels on the left HS
choiceDur = [0.25 1.8 0.25]; %duration of different choice windows in seconds (make sure to use durations that occur frequently enough)
stimDur = 1.5; %duration of stimulus sequence
sRate = 15; %sampling rate in Hz
groups = {'mSM', 'Fez', 'Plex' 'CSP'};
groupColor = {[1 0 0] [0 0 1] [0 0 0] [0 1 0]}; %colors for groups
animals = animals(1:10);

% this determines which regressors to use
dualRegs = {'firstAudStim' 'handleChoice' 'stimChoice' 'responseChoice' 'Grab', 'Lick'}; %choice, stimulus regressors
choiceRegs = {'handleChoice', 'stimChoice', 'responseChoice'};

% this determines the location of different areas
traceIdx = [85 390; 170 340; 250 155; 175 275]; %center of target areas
traceRad = [15, 15, 15, 15]; %radius of each area
traceLabels = {'Auditory' 'Parietal' 'Frontal' 'Somatosensory'};

% compute masks for each area, based on coordinates
nrTraces = size(traceIdx,1);
[xx,yy] = meshgrid(1:size(allenMask,2)/2,1:size(allenMask,1)); %isolate index for selected area
for x = 1 : nrTraces
    mask(:,:,x) = ~(hypot(xx - traceIdx(x,1), yy - traceIdx(x,2)) <= traceRad(x)); %circle masks
end

%% load beta maps
groupCnt = cell(1,length(groups));
dualBeta = cell(length(animals), length(dualRegs), 2);
diffBeta = cell(length(animals), length(dualRegs));
for iAnimals = 1 : length(animals)
    
    % check group ID
    for iGroups = 1 : length(groups)
        if contains(animals{iAnimals}, groups{iGroups})
            groupCnt{iGroups} = [groupCnt{iGroups} iAnimals];
        end
    end
    
    % load dual regressors
    for iRegs = 1 : length(dualRegs)
        load([cPath animals{iAnimals} '\blockData\' fileExt 'betaMaps_' trainingRange '_l' dualRegs{iRegs} '.mat'])
        dualBeta{iAnimals, iRegs, 1} = cBeta;
        dSize = size(cBeta,2);
        
        load([cPath animals{iAnimals} '\blockData\' fileExt 'betaMaps_' trainingRange '_r' dualRegs{iRegs} '.mat'])
        dualBeta{iAnimals, iRegs, 2} = cBeta;
        dSize = min([dSize, size(cBeta,2)]);
        
        % compute L/R difference
        diffBeta{iAnimals, iRegs} = dualBeta{iAnimals, iRegs, 1}(:,1:dSize) - dualBeta{iAnimals, iRegs, 2}(:,1:dSize); %keep L/R difference
    end
end

%% compute differential and single hemissphere results
[~,dSize] = cellfun(@size, diffBeta); %get sizes
dSize = min(dSize,[],1); %minimum regressor number
for x = 1 : size(diffBeta,1)
    for y = 1 : size(diffBeta,2)
        diffBeta{x,y} = diffBeta{x,y} (:,1:dSize(y));
    end
end

for iGroups = 1 : length(groups)
    for iRegs = 1 : length(dualRegs)
        meanDiff{iGroups, iRegs} = arrayShrink(nanmedian(cat(3,diffBeta{groupCnt{iGroups},iRegs}),3), allenMask, 'split');
        diffHS{iGroups, iRegs} = rateDisc_flipHS(meanDiff{iGroups, iRegs}, 'diff');
    end
end

%% forepaw figure
figure('Renderer','painters');
cIdx = strcmpi(dualRegs, 'Grab');
cRange = [0.00125 0.00125 0.00125 0.00125]; %color range
tRange =[7 25; 7 20; 7 20; 7 20];
for iGroups = 1 : length(groups)
    subplot(1,length(groups),iGroups);
    cData = smoothImg(nanmean(diffHS{iGroups,cIdx}(:,:,tRange(iGroups,1):tRange(iGroups,2)),3),1,5,5);
    imageScale(cData); axis image; colormap(colormap_blueblackred(256));
    caxis([-cRange(iGroups) cRange(iGroups)]);
    title([groups{iGroups} ' - Forepaw'])
    rateDisc_plotAllenOutline(gca,'L');
end

% make traces
figure('Renderer','painters');
cTraces = NaN(size(diffHS{1,cIdx},3), nrTraces, length(groups));
cRange = [100 60 60 60] .* 10^-4; %amplitude range
for iTraces = 1 : nrTraces
    for iGroups = 1 : length(groups)
        cTraces(:, iTraces, iGroups) = nanmean(arrayShrink(diffHS{iGroups,cIdx}, mask(:,:,iTraces), 'merge'),1);
    
        %plot traces
        subplot(nrTraces, length(groups), (iTraces-1).* length(groups) + iGroups)
        plot(cTraces(:,iTraces,iGroups), 'k', 'linewidth', 2); axis tight
        ylim([-4E-4 cRange(iGroups)]);  nhline(0, 'g'); nvline(7, 'g'); nvline(tRange(iGroups,:), 'color', [0.5 0.5 0.5]);
        title([groups{iGroups} ' - ' traceLabels{iTraces}]);
    end
end

%% stimulus figure
figure('Renderer','painters');
cIdx = strcmpi(dualRegs, 'firstAudStim');
cRange = [0.00025 0.0005 0.0005 0.0005]; %color range
tRange =[18 22; 18 22; 18 22; 18 22];
for iGroups = 1 : length(groups)
    subplot(1,length(groups),iGroups);
    cData = smoothImg(nanmean(diffHS{iGroups,cIdx}(:,:,tRange(iGroups,1):tRange(iGroups,2)),3),1,5,5);
    imageScale(cData); axis image; colormap(colormap_blueblackred(256));
    caxis([-cRange(iGroups) cRange(iGroups)]);
    title([groups{iGroups} ' - Stimulus'])
    rateDisc_plotAllenOutline(gca,'L');
end

% make traces
figure('Renderer','painters');
cTraces = NaN(size(diffHS{1,cIdx},3), nrTraces, length(groups));
cRange = [10 10 10 10] .* 10^-4; %amplitude range
for iTraces = 1 : nrTraces
    for iGroups = 1 : length(groups)
        cTraces(:, iTraces, iGroups) = nanmean(arrayShrink(diffHS{iGroups,cIdx}, mask(:,:,iTraces), 'merge'),1);
    
        %plot traces
        subplot(nrTraces, length(groups), (iTraces-1).* length(groups) + iGroups)
        plot(cTraces(:,iTraces,iGroups), 'k', 'linewidth', 2); axis tight
        ylim([-4E-4 cRange(iGroups)]);  nhline(0, 'g'); nvline(15, 'g'); nvline(tRange(iGroups,:), 'color', [0.5 0.5 0.5]);
        title([groups{iGroups} ' - ' traceLabels{iTraces}]);
    end
end

%% choice
choiceIdx = find(ismember(dualRegs, choiceRegs));
for iGroups = 1 : length(groups)
    for iRegs = 1 : length(choiceRegs)
        useIdx = squeeze(~isnan(nanmean(nanmean(diffHS{iGroups,choiceIdx(iRegs)},1),2)))'; %remove NaNs
        useIdx(find(useIdx,3,'last')) = false; %don't use last reqressors (can be very noisy)
        
        aChoice{iGroups,iRegs} = meanDiff{iGroups,choiceIdx(iRegs)}(:,:,useIdx);
        sHsChoice{iGroups,iRegs} = diffHS{iGroups,choiceIdx(iRegs)}(:,:,useIdx);
        choiceSegs(iGroups,iRegs) = sum(useIdx); %keep separation between task episodes
    end
end

%% show average for different choice episodes
choiceMovie = cell(1,length(groups));
temp = NaN(size(sHsChoice{1},1),size(sHsChoice{1},2), 1, 'single');
mChoice = NaN(size(sHsChoice{1},1),size(sHsChoice{1},2), 4, length(groups), 'single');
meanRange = 0.5; %duration of average in each segment (seconds)
for iGroups = 1 : length(groups)
    Cnt = 0;
    for iRegs = 1 : length(choiceRegs)
        Cnt = Cnt+1;
        cData = sHsChoice{iGroups,iRegs}(:,:,1:round(choiceDur(iRegs)*sRate));
        choiceMovie{iGroups} = cat(3,choiceMovie{iGroups},temp,cData); %make choice movie for traces
        
        % make average map
        if iRegs == 2 %stimulus regressor
            cIdx = 1 : min([size(cData,3) floor(meanRange*sRate)]);
            mChoice(:,:,Cnt,iGroups) = nanmean(cData(:,:,cIdx),3);
            Cnt = Cnt+1;
        end
        cIdx = max([1 size(cData,3)-floor(meanRange*sRate)]):size(cData,3);
        mChoice(:,:,Cnt,iGroups) = nanmean(cData(:,:,cIdx),3);
    end
end
    
% make figure
choiceLabels = {'Baseline' 'Stimulus' 'Delay' 'Response'};
nrRegs = size(mChoice,3);
cRange = [0.00075 0.00075 0.0015 0.0015];
for iGroups = 1 : size(mChoice,4)
    figure('Renderer','painters');
    for iRegs = 1 : nrRegs
        subplot(1, nrRegs, iRegs);
        imageScale(smoothImg(mChoice(:,:,iRegs,iGroups),1,3,3)); axis image
        colormap(colormap_blueblackred(256));
        caxis([-cRange(iGroups) cRange(iGroups)]); 
        title([groups{iGroups} ' - ' choiceLabels{iRegs}]);
    end
end
% choiceMovie = choiceMovie - mChoice(:,:,1); %subtract 'baseline'

%% make choice traces
figure('Renderer','painters');
cTraces = NaN(size(choiceMovie{1},3), nrTraces, length(groups));
cRange = [6 8 20 20] .* 10^-4; %amplitude range

for iTraces = 1 : nrTraces
    for iGroups = 1 : length(groups)
        cTraces(:, iTraces, iGroups) = nanmean(arrayShrink(choiceMovie{iGroups}, mask(:,:,iTraces), 'merge'),1);
    
        %plot traces
        subplot(nrTraces, length(groups), (iTraces-1).* length(groups) + iGroups)
        plot((cTraces(:,iTraces,iGroups)), 'k', 'linewidth', 2); axis square
        ylim([-6E-4 cRange(iGroups)]); xlim([1 size(choiceMovie{1},3)]);
        nhline(0, 'g');
        title([groups{iGroups} ' - ' traceLabels{iTraces}]);
    end
end


%% make overview figure
% figure('Renderer','painters');
% cRange = 0.001; %color range
% cGroup = 'CSP';
% cIDx = contains(groups,cGroup);
% for iRegs = 1 : length(choiceRegs)
%     subplot(1,length(choiceRegs),iRegs);
%     
%     a = smoothImg(mChoice{cIDx}(:,:,iRegs),1,5,5);
%     imageScale(a); axis image; colormap(colormap_blueblackred(256));
%     caxis([-cRange cRange]);     
%     title([cGroup ' - ' choiceRegs{iRegs}]);
%     niceFigure(gca);
% end
% 
% %make handle vs stim figure
% figure('Renderer','painters');
% cRange = 0.00025; %color range
% ind1 = ismember(choiceRegs,'handleChoice');
% ind2 = ismember(choiceRegs,'stimChoice');
% 
% subplot(1,2,1);
% cData = smoothImg(mChoice(:,:,ind2),1,5,5) - smoothImg(mChoice(:,:,ind1),1,5,5);
% imageScale(cData); axis image; colormap(colormap_blueblackred(256));
% caxis([-cRange cRange]); 
% 
% % 
% figure('Renderer','painters');
% cRange = 0.0002; %color range
% subplot(1,2,1);
% stimStart = round(choiceDur(1)*sRate)+10; %start of stim epoch
% 
% cData = nanmean(choiceMovie(:,:,stimStart:stimStart+5),3);
% imageScale(smoothImg(cData,1,5,5)); caxis([-cRange cRange]); axis image; colormap(colormap_blueblackred(256));


%% stimulus
% figure('Renderer','painters');
% for iGroups = 1 : length(groups)
%     subplot(1,length(groups),iGroups);
%     cRange = 0.00025; %color range
%     cData = smoothImg(nanmean(sHsStim{iGroups}(:,:,18:21),3),1,5,5);
%     imageScale(cData); axis image; colormap(colormap_blueblackred(256));
%     caxis([-cRange cRange]);
%     title('First audio pulse')
%     rateDisc_plotAllenOutline(gca,'L');
% end

% subplot(1,2,2);
% cRange = 0.0001; %color range
% cData = smoothImg(nanmean(sHsStim{2}(:,:,1:2),3),1,5,5);
% imageScale(cData); caxis([-cRange cRange]); axis image; colormap(colormap_blueblackred(256));
% title('Other audio pulses')

% 
% 
% %%
% lChoice1 = smoothCol(lChoice,2);
% rChoice1 = smoothCol(rChoice,2);
% 
% lhandleChoice1 = smoothCol(lhandleChoice,2);
% rhandleChoice1 = smoothCol(rhandleChoice,2);
% 
% lStim1 = smoothCol(lStim,2);
% rStim1 = smoothCol(rStim,2);
% 
% aChoice = arrayShrink(nanmean(lChoice - rChoice,3), allenMask, 'split');
% ahandleChoice = arrayShrink(nanmean(lhandleChoice - rhandleChoice,3), allenMask, 'split');
% aStim = arrayShrink(nanmean(lStim - rStim,3), allenMask, 'split');
