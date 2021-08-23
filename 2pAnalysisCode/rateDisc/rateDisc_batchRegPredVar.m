% function rateDisc_batchRegPredVar
% code to run over all recordings collect predicted variance for individual
% regressors and regressor groups.

cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
savePath = 'R:\BpodImager\predVar\'; %local save path
disp(cPath);

load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
% groups = {'mSM' 'Fez' 'Plex'};
groups = {'mSM' 'Fez' 'Plex' 'CSP'};
reload = false; %reload data from server

% get animal info
[dataOverview, motorLabels, sensorLabels, cogLabels, ~, ~, ~, ~, trainDates, allRegs] = rateDiscRecordings;
trainingRange = 'audioDisc'; %use this to run analysis only in a certain range of training
animals = dataOverview(:,1);
animals = animals(1:10);
fileExt = 'org';

%% load data
if reload
    % get data
    allVar = NaN(sum(~allenMask(:)), length(allRegs), length(animals), 2, 'single');
    fullVar = NaN(sum(~allenMask(:)), length(animals), 5, 'single');
    for iAnimals = 1 : length(animals)
        
        % get recs
        cAnimal = animals{iAnimals};
        recs = dir([cPath cAnimal filesep 'SpatialDisc' filesep]);
        recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs, trainingRange); %this sorts recordings by date
        
        cData = NaN(sum(~allenMask(:)), length(allRegs), 2, length(recs), 'single');
        cMaps = NaN(sum(~allenMask(:)), 3, length(recs), 'single');
        for iRecs = 1 : length(recs)
            
            try
                disp([cAnimal ' - ' recs(iRecs).name]);
                fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %session to predicted variance results
                
                load([fPath 'opts2.mat'],'opts')
                load([fPath fileExt 'fullcorr.mat'],'fullMap')
                load([fPath 'taskregData.mat'],'taskMap')
                load([fPath 'motorregData.mat'],'motorMap')
                
                try
                    load([fPath 'mask.mat'],'mask')
                catch
                    matFile = ([fPath 'rsVc.mat']);
                    mask = isnan(matFile.U(:,:,1));
                    save([fPath 'mask.mat'],'mask')
                end
                
                fullMap = arrayShrink(fullMap, mask, 'split');
                fullMap = alignAllenTransIm(single(fullMap.^2),opts.transParams);
                cMaps(:, 1, iRecs) = arrayShrink(fullMap, allenMask, 'merge');
                taskMap = arrayShrink(taskMap, mask, 'split');
                taskMap = alignAllenTransIm(single(taskMap.^2),opts.transParams);
                cMaps(:, 2, iRecs) = arrayShrink(taskMap, allenMask, 'merge');
                motorMap = arrayShrink(motorMap, mask, 'split');
                motorMap = alignAllenTransIm(single(motorMap.^2),opts.transParams);
                cMaps(:, 3, iRecs) = arrayShrink(motorMap, allenMask, 'merge');
                
                fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep 'predVariance' filesep];
                for iRegs = 1 : length(allRegs)
                    
                    cLabel = [fileExt allRegs{iRegs}];
                    load([fPath 'shCurrent' filesep cLabel 'corr.mat'], 'cMap');
                    cData(:, iRegs, 1, iRecs) = cMap.^2;
                    
                    if ~strcmpi(cLabel, [fileExt 'full'])
                        load([fPath 'shOther' filesep cLabel 'corr.mat'], 'cMap');
                        cData(:, iRegs, 2, iRecs) = cMap.^2;
                    end
                    
                    if ~exist('oMotorLabels','var')
                        load([fPath fileExt 'oMotorLabels.mat'], 'oMotorLabels');
                    end
                    
                end
            catch
                cData(:,:,:,iRecs) = NaN(size(cData,1),size(cData,2),size(cData,3));
                cMaps(:,:,iRecs) = NaN(size(cMaps,1),size(cMaps,2));
                fprintf('Couldnt load %s\n', [fPath 'shCurrent' filesep cLabel 'corr.mat']);
            end
        end
        allVar(:, :, iAnimals, :) = nanmean(cData,4); clear cData
        fullVar(:, iAnimals, 1:3) = nanmean(cMaps,3);
        fullVar(:, iAnimals, 4) = nanmean(cMaps(:,1,:) - cMaps(:,2,:),3); %unique movement
        fullVar(:, iAnimals, 5) = nanmean(cMaps(:,1,:) - cMaps(:,3,:),3); %unique task
    end
    % save data
    save([savePath 'regVariance.mat'], 'allVar', 'fullVar','oMotorLabels','-v7.3')
else
    load([savePath 'regVariance.mat'], 'allVar', 'fullVar','oMotorLabels')
end

%% collect data for groups
selVar = cell(1, length(groups));
selMaps = cell(1, length(groups));
for x = 1 : length(groups)
    cIdx = find(contains(animals,groups{x}))'; %animals that belong to current group
    selVar{x} = allVar(:,:,cIdx,:); %data for current group
    selMaps{x} = fullVar(:,cIdx,:); %data for current group
end

%% make some figures
idx = zeros(1,length(allRegs)-1);
idx(1:10) = 1;
idx(11:12) = 2;
idx(13:end) = 3;
idxGroup{1} = 1:sum(idx == 1);
idxGroup{2} = sum(idx == 1)+1:sum(ismember(idx, 1:2));
idxGroup{3} = sum(ismember(idx, 1:2))+1:sum(ismember(idx, 1:3));

% show explained variance for each variable
for x = 1 : length(groups)
    figure('Renderer','painters');
    cData = squeeze(nanmean(selVar{x},1));
    [ax, ~] = regressorBoxPlot(cData(2:end,:,2)', allRegs(2:end), 5, subplot(2,1,1), [0 1 0], idxGroup, [0 0.4]); title(groups{x});
    
    regressorBoxPlot(bsxfun(@minus,cData(1,:,1),cData(2:end,:,1))', allRegs(2:end), 5, subplot(2,1,2), [25 111 61]/255, idxGroup,[0 0.05]); title(groups{x});
end

%% show variable matrix across mice
noUseLabels = {'full' 'lfirstAudStim' 'rfirstAudStim' 'lAudStim' 'rAudStim' 'handleChoice' ...
    'stimChoice' 'respChoice' 'handleSound' 'video'};

idx = zeros(1,length(allRegs));
idx(~ismember(allRegs,[motorLabels, oMotorLabels])) = 1; %task
idx(ismember(allRegs,[motorLabels, oMotorLabels])) = 2; %movement
idx(ismember(allRegs,noUseLabels)) = 0; %dont use these
vidReg = ismember(allRegs,'video');

figure('Renderer','painters');
subplot(1,2,1)
a = [];
for x = 1 : length(groups)
    b = squeeze(nanmean(selVar{x}(:,idx==1,:,2),1));
    %     b = bsxfun(@rdivide,b,squeeze(nanmean(selVar{x}(:,1,:,1),1))');
%     b = bsxfun(@minus,squeeze(nanmean(selVar{x}(:,1,:,1),1))',b);
%         b = bsxfun(@rdivide,b,squeeze(nanmean(selVar{x}(:,vidReg,:,1),1))'); %divide non-video model
    a = [a b];
end
[u,i] = sort(nanmedian(a,2), 'descend');
ax = imagesc(a(i,:)'); title('task - cvR^2');
ax.Parent.XTick = 1:sum(idx==1);
cLabels = allRegs(ismember(idx,1));
ax.Parent.XTickLabels = cLabels(i);
ax.Parent.XTickLabelRotation = 45;
ax.Parent.YTickLabels = animals;
axis image;
caxis([0 0.2]);
ax.Parent.YTick = 1:length(animals);
ax.Parent.YTickLabels = animals;
ax.Parent.TickLength = [0 0];
niceFigure(ax.Parent)

subplot(1,2,2)
a = [];
for x = 1 : length(groups)
    b = squeeze(nanmean(selVar{x}(:,ismember(idx,2:3),:,2),1));
    %     b = bsxfun(@rdivide,b,squeeze(nanmean(selVar{x}(:,1,:,1),1))');
%     b = bsxfun(@minus,squeeze(nanmean(selVar{x}(:,1,:,1),1))',b);
%         b = bsxfun(@rdivide,b,squeeze(nanmean(selVar{x}(:,vidReg,:,1),1))'); %divide non-video model
    a = [a b];
end
[u,i] = sort(nanmedian(a,2),'descend');
ax = imagesc(a(i,:)'); title('movement - cvR^2');
ax.Parent.XTick = 1:sum(ismember(idx,[2 3]));
cLabels = allRegs(ismember(idx,[2 3]));
ax.Parent.XTickLabels = cLabels(i);
ax.Parent.XTickLabelRotation = 45;
ax.Parent.YTickLabels = animals;
axis image;
caxis([0 0.2]);
ax.Parent.YTick = 1:length(animals);
ax.Parent.YTickLabels = animals;
ax.Parent.TickLength = [0 0];
niceFigure(ax.Parent)
colormap(viridis(256));

%% explained variance maps
figure('Renderer','painters');
cRange = [0.7 0.45 0.06 0.01];
for x = 1 : length(groups)
    
    subplot(length(groups),4,(x-1)*4+1)
    cData = nanmean(selMaps{x}(:,:,1),2); %full model
    imageScale(arrayShrink(cData, allenMask, 'split'));
    axis image; colormap(inferno(256)); caxis([0 cRange(1)]);
    title('Full model cVR^2');
    
    subplot(length(groups),4,(x-1)*4+2)
    cData = nanmean(selMaps{x}(:,:,4),2); %move model
    imageScale(arrayShrink(cData, allenMask, 'split'));
    axis image; colormap(inferno(256)); caxis([0 cRange(2)]);
    title('Movement model dVR^2');
    
    subplot(length(groups),4,(x-1)*4+3)
    cData = nanmean(selMaps{x}(:,:,1),2) - nanmean(selMaps{x}(:,:,4),2) - nanmean(selMaps{x}(:,:,5),2); %shared move/task model
    imageScale(arrayShrink(cData, allenMask, 'split'));
    axis image; colormap(inferno(256)); caxis([0 cRange(3)]);
    title('Shared task/movement dVR^2');
    
    subplot(length(groups),4,(x-1)*4+4)
    cData = nanmean(selMaps{x}(:,:,5),2); %task model
    imageScale(arrayShrink(cData, allenMask, 'split'));
    axis image; colormap(inferno(256)); caxis([0 cRange(4)]);
    title('Task model dVR^2');
    
end

%% single variable maps
Cnt = 0;
for x = 1 : length(groups)
    Cnt = max([Cnt, size(selVar{x},3)]); %find largest group
end
Cnt = Cnt + 1; %add average mouse for each group
clear h

% for iRegs = 1 : length(cRegs)
h = figure('Renderer','painters');
cReg = 'audioStim';
for x = 1 : length(groups)
    
    for y = 1 : size(selVar{x},3)
        cIdx = strcmpi(allRegs,cReg); %current regressors
        cMap = squeeze(selVar{x}(:,1,y,1)) - squeeze(selVar{x}(:,cIdx,y,1)); %subtract 'shCurrent' map for current regressor from full model
        
        subplot(length(groups),Cnt,(x-1)*Cnt + y);
        [~, cRange] = imageScale(arrayShrink(cMap,allenMask)); title([groups{x} num2str(y) ' - ' cReg]);
        caxis([0 cRange(2)]);
        colormap(inferno(256));
    end
    
    %average mouse
    cMap = nanmean(squeeze(selVar{x}(:,1,:,1)) - squeeze(selVar{x}(:,cIdx,:,1)), 2); %subtract 'shCurrent' map for current regressor from full model
    subplot(length(groups),Cnt,(x-1)*Cnt + y + 1);
    [~, cRange] = imageScale(arrayShrink(cMap,allenMask)); title([groups{x} 'avg - ' cReg]);
    colormap(inferno(256)); drawnow;
    caxis([0 cRange(2)]);
end


%%
% 
% % bar plots explained variance
% for x = 1
%     
%     % stack unique task variance (line 1), shared varance between task and motor (line 2) and unique motor variance (line 3) together
%     cData = [nanmean(selMaps{x}(:,:,5)); nanmean(selMaps{x}(:,:,3))-nanmean(selMaps{x}(:,:,4)); nanmean(selMaps{x}(:,:,4))];
%     
%     figure('Renderer','painters');
%     pie(flipud(double(nanmean(cData,2)./mean(nanmean(selMaps{x}(:,:,1))))))
%     %     pie(flipud(double(nanmean(cData,2))))
%     
% end

