% code to analyze locaNMF results and check changes over time
cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
tPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\BpodImager\Animals\'; %path to churchlandNAS
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
% animals = animals([1:7 9:10]);
animals = animals(1);
load('allenDorsalMapSM.mat')
[xRange, yRange] = rateDisc_maskRange(dorsalMaps.allenMask); % get inner range of allenMask
allenMask = dorsalMaps.allenMask(yRange,xRange);
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
dimCnt = 15; %number of dimensions from each area that are considered for dimension correlations
groups = {'mSM' 'Fez' 'Plex' 'CSP'};

%% get regions
load('trimMap.mat','trimMap')
trimMap = trimMap(yRange,xRange);
useAreas = unique(arrayShrink(trimMap,allenMask)); %areas that are in the filed of view

% define larger regions
regLabels = {'ACC' 'M2' 'M1' 'SSf' 'SSb' 'Aud' 'PPC' 'V1' 'V2' 'RS'}; %labels for different regions
% reg{1} = find(ismember(dorsalMaps.labelsSplit,'MOB')); %olfactory
reg{1} = find(contains(dorsalMaps.labelsSplit,'ACAd')); %ACC
reg{2} = find(contains(dorsalMaps.labelsSplit,'MOs')); %secondary motor
reg{3} = find(contains(dorsalMaps.labelsSplit,'MOp')); %primary motor
reg{4} = find(contains(dorsalMaps.labelsSplit,'SSp-bfd') | contains(dorsalMaps.labelsSplit,'SSp-m')); %somatosensory face (whisker and mouth)
reg{5} = find(contains(dorsalMaps.labelsSplit,'SS')); %somatosensory body
reg{5}(ismember(reg{5},reg{4})) = []; % exclude SSface
reg{6} = find(contains(dorsalMaps.labelsSplit,'AUD')); %auditory
reg{7} = find(contains(dorsalMaps.labelsSplit,'VISal') | contains(dorsalMaps.labelsSplit,'VISrl')); %parietal (using AL and RL)
reg{8} = find(strcmpi(dorsalMaps.labelsSplit,'VISp')); %primary visual
reg{9} = find(contains(dorsalMaps.labelsSplit,'VIS')); %secondary visual
reg{9}(ismember(reg{9},reg{8}) | ismember(reg{9},reg{7})) = []; % exclude V1 and parietal from V2
reg{10} = find(contains(dorsalMaps.labelsSplit,'RSPagl')); %retrosplinal cortex

allDims = cat(1,reg{:}); %all area IDs that will be used
allDims = allDims(ismember(allDims, useAreas)); %only use areas within field of view

%% run over animals
dimCorr = NaN(length(allDims) * dimCnt, length(allDims) * dimCnt, length(animals), 'single'); %sorted correlation map between dimensions
regCorr = NaN(length(reg), length(reg), length(animals), 'single'); %sorted correlation map between dimensions
allA = cell(1,length(animals));
for iAnimals = 1 : length(animals)
    
    cAnimal = animals{iAnimals}; % current animal
    fprintf('Current animal: %s\n', cAnimal);
    bPath = [cPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    tbPath = [tPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    if ~exist(tbPath,'dir')
        mkdir(tbPath);
    end
    load([tbPath 'mask_' trainingRange '.mat'], 'mask');
    load([tbPath 'locanmf_' cAnimal '_loc50_' trainingRange '.mat'], 'A', 'C', 'areas');
    allA{iAnimals} = arrayShrink(A,allenMask,'merge');
    dimNr(iAnimals) = size(A,3);
    [~,b] = max(diff(areas));
    areas(b+1:end) = 256 - areas(b+1:end); %match areas with dorsalmap IDs
    
    %make region correlation map
    [regCorr(:,:,iAnimals), dimCorr(:,:,iAnimals), B, regIdx] = rateDisc_regionCorr(A, C, areas, reg, allDims, dimCnt, useAreas, dorsalMaps);
    
    %compute explained variance for each dimension
    A = arrayShrink(A,mask,'merge');
    covC = diag(cov(bsxfun(@minus, C, mean(C,2))'));  % S x 1

    varC = NaN(size(A), 'single');
    for x = 1 : size(A,1)
        varC(x,:) = A(x,:) .* covC' .* A(x,:);
    end
    
    C = reshape(C,size(C,1),[]);
    varInDimD = NaN(1, size(A, 2));
    for d = 1:size(C, 1)
      covV1 = cov(C(:,d)');  % S x S
      varP1 = sum((A * covV1) .* A, 2)';  % 1 x P
      varInDimD(d) = sum(varP1);
    end
end

%% make regional prediction figure
figure
for iGroups = 1 : length(groups)
    cGroup = (contains(animals', groups{iGroups}));
    subplot(2,2,iGroups);
    ax = imagesc(nanmean(real(regCorr(:,:,cGroup)),3)); 
    title(['Predicted variance - R^2: ' groups{iGroups}]);
    ax.Parent.XTick = 1:size(regCorr,2);
    ax.Parent.XTickLabelRotation = 45;
    ax.Parent.YTickLabels = regLabels;
    ax.Parent.YTick = 1:size(regCorr,2);
    ax.Parent.TickLength = [0 0];
    ylabel('Predictor regions');
    niceFigure(ax.Parent)
    axis square; colormap('inferno');
    caxis([0.1 0.9]);
    hline([2.5 6.5], '--w') %outline for motor/ss regions
    vline([2.5 6.5], '--w') %outline for motor/ss regions
    hline([7.5 10.5], '--g') %outline for vision regions
    vline([7.5 10.5], '--g') %outline for vision regions
    niceFigure(ax.Parent)
    if iGroups > 2
        ax.Parent.XTickLabels = regLabels;
        xlabel('Predicted regions');
    else
        ax.Parent.XTickLabels = [];
    end
end

%% make component correlation maps
lineIdx = repmat(allDims',dimCnt,1);
lineIdx = reshape(lineIdx, 1, []);
a = mean(dimCorr,3);
rejIdx = isnan(nanmean(a,2));
lineIdx(rejIdx) = [];
a(rejIdx,:) = [];
a(:,rejIdx) = [];

clear regMean regOff
for iRegs = 1 : length(reg)
    if any(ismember(lineIdx, reg{iRegs}))
        regMean(iRegs) = nanmean(find(ismember(lineIdx, reg{iRegs})));
        regOff(iRegs) = find(ismember(lineIdx, reg{iRegs}),1, 'last');
    else
        regMean(iRegs) = NaN;
    end
end

figure
b = tril(smoothImg(nanmean(real((a)),3),2,1.1));
b(b == 0) = NaN;
ax = imagesc(b); title('Predicted variance - R^2');
ax.Parent.XTick = regMean;
ax.Parent.XTickLabelRotation = 45;
ax.Parent.XTickLabels = regLabels(~isnan(regMean));
ax.Parent.YTickLabels = regLabels(~isnan(regMean));
ax.Parent.YTick = regMean;
ax.Parent.TickLength = [0 0];
set(ax,'AlphaData',~isnan(ax.CData)); %make NaNs transparent.
xlabel('Predicted regions');
ylabel('Predictor regions');
axis square; colormap('parula');
caxis([-0.03 0.03]);
hline(regOff(1:end-1)+0.5, '--w');
vline(regOff(1:end-1)+0.5, '--w');
niceFigure(ax.Parent)

%% check similarity of spatial components
cGroup = find(contains(animals', groups{1}));
a = cat(2,allA{:});
b = a(~isnan(mean(a,2)),:);
b = corrcoef(b(:,nanmean(b,1) ~= 0));

%%
figure;
ax = imagesc(b(:,1:sum(dimNr(1)))'); colorbar; axis image;
% ax = imagesc(b(:,sum(dimNr(1)):sum(dimNr(1:2)))'); colorbar; axis image;
title('Correlation across animals');
ylabel(['Spatial components - ' animals{cGroup(1)}]);
xlabel(['Spatial components - ' groups{1}]);
niceFigure(ax.Parent);

%%



