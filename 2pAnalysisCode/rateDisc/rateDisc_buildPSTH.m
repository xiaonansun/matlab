function rateDisc_buildPSTH(group,trainingRange)

if ispc
    %     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
    %     cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
    cPath = 'Q:\BpodImager\Animals\'; %local data path
else
    %     cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/'; %path to grid
end

if ~exist('trainingRange','var') || isempty(trainingRange)
    trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
end

savePath = 'R:\BpodImager\trialAvg\'; %local save path
load('allenDorsalMapSM.mat', 'dorsalMaps')
allenMask = dorsalMaps.allenMask;
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
animals = animals(contains(animals,group));
reload = false;
opts.frameRate = 15;
opts.preStim = 3;
segIdx = [1 0.75 1.25 0.75 0.75];
segFrames = cumsum(floor(segIdx * opts.frameRate)); %max nr of frames per segment
frames = 90; %nr of frames per trial
baseWin = floor(opts.frameRate/2); %nr of baseline frames that are set to 0
dWin = floor(opts.frameRate/2); %nr of frames in stimulus/delay that are used to compute d'

% this determines the location of different areas
traceIdx = [85 390; 170 340; 250 155]; %center of target areas (left HS)
% traceIdx(:,1) = 586-traceIdx(:,1); %switch to right HS
traceRad = [15, 15, 15]; %radius of each area
traceLabels = {'Sensory' 'Parietal' 'Frontal'};

% compute masks for each area, based on coordinates
nrTraces = size(traceIdx,1);
[xx,yy] = meshgrid(1:size(allenMask,2),1:size(allenMask,1)); %isolate index for selected area
for x = 1 : nrTraces
    mask(:,:,x) = ~(hypot(xx - traceIdx(x,1), yy - traceIdx(x,2)) <= traceRad(x)); %circle masks
end

%% load data
Cnt = 0;
figure('Renderer','painters');
allPSTH = NaN(sum(~allenMask(:)), segFrames(end), 5, length(animals), 'single');
allPSTHvar = NaN(sum(~allenMask(:)), 2, length(animals), 'single');
allPrime = NaN(sum(~allenMask(:)), length(segFrames), 2, length(animals), 'single');
for iAnimals = 1:length(animals)
    clear cPSTH cVar dPrime
    cReload = reload;
    if ~reload
        
        try
            % load psth data
            load([savePath 'allPSTH_' trainingRange '_' animals{iAnimals}], 'cPSTH', 'cVar', 'dPrime')
        catch
            disp(['Couldnt find data for mouse ' animals{iAnimals} '. Regenerate from raw data instead.'])
            cReload = true;
        end
    end
    
    if cReload
        disp(animals{iAnimals});
        recs = rateDisc_getRecsForAnimal(animals{iAnimals},trainingRange, cPath);% end
        cPSTH = NaN(sum(~allenMask(:)), segFrames(end), 5, length(recs), 'single');
        cVar = NaN(sum(~allenMask(:)), 2, length(recs), 'single');
        dPrime = NaN(sum(~allenMask(:)), length(segFrames), 2, length(recs), 'single');
        for iRecs = 1 : length(recs)
            try
                % get imaging data
                fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep];                
                try
                    load([fPath 'rsVc.mat'], 'Vc', 'bTrials');
                catch
                    load([fPath 'Vc.mat'], 'Vc', 'bTrials');
                    if size(Vc,2) > 105
                        Vc = [];
                    end
                end
                Vc = Vc(:, 1:frames, :);
                load([fPath 'alignU.mat'], 'U');
                U = rateDisc_removeOutline(U, 10);
                U = arrayShrink(U,allenMask,'merge');
                
                % get behavior and different trial averages
                bhvFile = dir([fPath animals{iAnimals} '_SpatialDisc*.mat']);
                load([fPath bhvFile(1).name], 'SessionData');
                bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset
                choiceIdx = bhv.ResponseSide == 1; %left response side
                stimIdx = bhv.CorrectSide == 1; %left side is correct
                singleIdx = bhv.DistStim == 0; %detection trials
                
                % produce variance maps
                a = reshape(Vc, size(Vc,1), []); %all data
                covV = cov(a(:,~isnan(nanmean(a,1)))');
                cVar(:,1,iRecs) = dot((U*covV)', U'); % 1 x P
                
                covV = cov(nanmean(Vc,3)'); %variance for trial average
                cVar(:,2,iRecs) = dot((U*covV)', U'); % 1 x P
                
                % get PSTH movies
                Vc = rateDisc_getBhvRealignment(Vc, bhv, segFrames, opts); %align to different trial segments
                baseAvg = nanmean(nanmean(Vc(:,1:baseWin,:),2),3);

                Vc = Vc - baseAvg; %remove baseline for PSTH
                cPSTH(:,:,1,iRecs) = U * nanmean(Vc,3); %all trial PSTH
                cPSTH(:,:,2,iRecs) = U * nanmean(Vc(:,:,choiceIdx),3); %left choice PSTH
                cPSTH(:,:,3,iRecs) = U * nanmean(Vc(:,:,~choiceIdx),3); %right choice PSTH
                cPSTH(:,:,4,iRecs) = U * nanmean(Vc(:,:,stimIdx & singleIdx),3); %left single stim PSTH
                cPSTH(:,:,5,iRecs) = U * nanmean(Vc(:,:,~stimIdx & singleIdx),3); %right single stim PSTH
                
                % compute variance maps for stimulus and delay episodes
                cSegs = [0 segFrames];
                for iSegs = 1 : length(cSegs)-1
                    
                    % choice (left vs right)
                    cData = squeeze(nanmean(Vc(:,cSegs(iSegs)+1:cSegs(iSegs)+dWin,choiceIdx),2));
                    cData = cData(:, ~isnan(cData(1,:)));
                    cVar1 = dot((U*cov(cData'))', U');
                    cMean1 = nanmean(cData,2);
                    cData = squeeze(nanmean(Vc(:,cSegs(iSegs)+1:cSegs(iSegs)+dWin,~choiceIdx),2));
                    cData = cData(:, ~isnan(cData(1,:)));
                    cVar2 = dot((U*cov(cData'))', U');
                    cMean2 = nanmean(cData,2);
                    
                    dPrime(:,iSegs,1,iRecs) = (U * (cMean1 - cMean2)) ./ sqrt((cVar1+cVar2)/2)'; %choice difference / standard deviations
                    
                    % stim (left vs right)
                    cData = squeeze(nanmean(Vc(:,cSegs(iSegs)+1:cSegs(iSegs)+dWin,stimIdx & singleIdx),2));
                    cData = cData(:, ~isnan(cData(1,:)));
                    cVar1 = dot((U*cov(cData'))', U');
                    cMean1 = nanmean(cData,2);
                    cData = squeeze(nanmean(Vc(:,cSegs(iSegs)+1:cSegs(iSegs)+dWin,~stimIdx & singleIdx),2));
                    cData = cData(:, ~isnan(cData(1,:)));
                    cVar2 = dot((U*cov(cData'))', U');
                    cMean2 = nanmean(cData,2);
                    
                    dPrime(:,iSegs,2,iRecs) = (U * (cMean1 - cMean2)) ./ sqrt((cVar1+cVar2)/2)'; %choice difference / standard deviations
                end
                fprintf('Done: %s - %s !!\n', animals{iAnimals}, recs(iRecs).name)
                
            catch ME
                fprintf('!! Warning: Failed to run model %s - %s !!\n', animals{iAnimals}, recs(iRecs).name)
                disp(ME.message);
            end
        end
        
        %% save result
        if ~exist(savePath,'dir')
            mkdir(savePath);
        end
        save([savePath 'allPSTH_' trainingRange '_' animals{iAnimals}], 'cPSTH', 'cVar', 'dPrime', 'segFrames', 'segIdx', '-v7.3');
    end
    
    % keep data. Looping saves some memory here.
    for iSegs = 1 : length(segFrames)
        for iFrames = 1 : size(cPSTH,2)
            allPSTH(:, iFrames, iSegs, iAnimals) = nanmedian(cPSTH(:,iFrames,iSegs,:),4);
        end
        allPrime(:, iSegs, :, iAnimals) = nanmedian(dPrime(:,iSegs,:,:),4);
    end
    allPSTHvar(:, :, iAnimals) = nanmedian(cVar,3);
    
    %% show variance figure
    subplot(length(animals),2,Cnt + 1);
    imageScale(arrayShrink(allPSTHvar(:,1,iAnimals), allenMask, 'split'));
    title([animals{iAnimals} ' - allVar']); drawnow;
    subplot(length(animals),2,Cnt + 2);
    imageScale(arrayShrink(allPSTHvar(:,2,iAnimals), allenMask, 'split'));
    colormap(inferno(256));
    title([animals{iAnimals} ' - trialVar']); drawnow;
    
    Cnt = Cnt + 2;
end

%% show variance figure
figure('Renderer','painters');
subplot(1,2,1);
cImg = imageScale(arrayShrink(nanmean(allPSTHvar(:,1,:),3), allenMask, 'split'));
title([group ' - allVar']); cImg.Parent.CLim(1) = 0;

subplot(1,2,2);
cImg = imageScale(arrayShrink(nanmean(allPSTHvar(:,2,:),3), allenMask, 'split'));
colormap(inferno(256)); cImg.Parent.CLim(1) = 0;
title([group ' - trialVar']);

%% compute PSTHs
figure('Renderer','painters');
% cRange = 0.0075; %color range
cRange = 0.02; %color range
cIdx = [0 segFrames];
segLabels = {'Baseline' 'Initiate' 'Stimulus' 'Delay' 'Response'};
cType = 1; %PSTH type
for iSegs = 1 : length(segFrames)
    subplot(2, 3, iSegs);
%     subplot(length(segFrames),1, iSegs);
    
%     cData1 = arrayShrink(nanmean(nanmean(allPSTH(:,cIdx(iSegs)+1 : cIdx(iSegs+1), 2, :), 4), 2), allenMask, 'split');
%     cData2 = arrayShrink(nanmean(nanmean(allPSTH(:,cIdx(iSegs)+1 : cIdx(iSegs+1), 3, :), 4), 2), allenMask, 'split');
%     cData = arrayShrink(nanmean(nanmean(allPSTH(:,cIdx(iSegs)+1 : cIdx(iSegs+1), cType, :), 4), 2), allenMask, 'split');
    cData = arrayShrink(nanmean(nanmean(allPSTH(:,cIdx(iSegs)+1 : cIdx(iSegs)+dWin, cType, :), 4), 2), allenMask, 'split');

    imageScale(smoothImg(cData)); axis image; colormap(colormap_blueblackred(256));
    caxis([-cRange cRange]);
    title([group ' - ' segLabels{iSegs}])
    rateDisc_plotAllenOutline(gca, 'R');
end
hold; plot(586-traceIdx(:,1),traceIdx(:,2), 'ko', 'linewidth',3); %show trace locations

%% make PSTH traces
figure('Renderer','painters');
yRange = 0.04; %y-axis range
cData = arrayShrink(squeeze(allPSTH(:,:, cType, :)), allenMask, 'split');
for iTraces = 1 : nrTraces
    
    cTraces = squeeze(nanmean(arrayShrink(cData, mask(:,:,iTraces), 'merge'),1));
    cTraces = insertColumns(cTraces', segFrames(2:end)+1, NaN)';
%     cTraces(1:9,:) = [];
    
    %plot traces
    subplot(nrTraces, 1, iTraces)
    plot(cTraces, 'linewidth', 1); hold on;  axis square
    plot(nanmean(cTraces,2), 'color','k', 'linewidth', 2);
    title([group ' - ' traceLabels{iTraces}]);
    ylim([-yRange/2 yRange]);
    nvline([segFrames(1); find(isnan(cTraces(:,1)))]+1, 'r', 'linewidth', 1); nhline(0,'k');
end

%% show dprime map
iArea = [3, 4];
segLabel = {'Stim' 'Delay'};
figure('Renderer','painters');
cRange = 0.2; %color range
Cnt = 0;
for x = 1:2
    subplot(1,4,Cnt+1);
    cData = arrayShrink(nanmean(allPrime(:,iArea(x),1,:),4), allenMask, 'split');
    cData = rateDisc_flipHS(cData,'diff');
    cImg = imageScale(smoothImg(cData, 1, 3, 3),cRange);
    title([group ' - d": Choice - ' segLabel{x}]);
    
    subplot(1,4,Cnt+2);
    cData = arrayShrink(nanmean(allPrime(:,iArea(x),2,:),4), allenMask, 'split');
    cData = rateDisc_flipHS(cData,'diff');
    cImg = imageScale(smoothImg(cData, 1, 3, 3),cRange);
    colormap(viridis(256));
    title([group ' - d": Sensory - ' segLabel{x}]);
    
    Cnt = Cnt+2;
end

%% make L/R difference traces
figure('Renderer','painters');
yRange = 0.0035; %y-axis range
cIdx = {[2 3], [4 5]}; %choice on first run, stimulus on second
cColor = {'k', 'g'};
for y = 1 : 2
    cData1 = arrayShrink(squeeze(allPSTH(:,:, cIdx{y}(1), :)), allenMask, 'split');
    cData1 = (rateDisc_flipHS(cData1,'diff'))./2;
    cData2 = arrayShrink(squeeze(allPSTH(:,:, cIdx{y}(2), :)), allenMask, 'split');
    cData2 = (rateDisc_flipHS(cData2,'diff'))./2;
    for iTraces = 1 : nrTraces
        
        cTraces = squeeze(nanmean(arrayShrink(cData1 - cData2, mask(1:size(cData,1),1:size(cData,2),iTraces), 'merge'),1));
        cTraces = insertColumns(cTraces', segFrames(2:end)+1, NaN)';
        
        %plot traces
%         subplot(2, nrTraces, (y-1)*nrTraces + iTraces)
        subplot(nrTraces, 1, iTraces)
        plot(nanmean(cTraces,2), 'color', cColor{y}, 'linewidth', 2); hold on;  axis square
%         stdshade(cTraces',[],cColor{y},0.5,1); hold on;  axis square
        title([group ' - ' traceLabels{iTraces}]);
        ylim([-yRange yRange]); xlim([segFrames(1) segFrames(end-1)]+1);
        nvline([segFrames(1); find(isnan(cTraces(:,1)))]+1, 'r', 'linewidth', 1); nhline(0,'k');
    end
end

