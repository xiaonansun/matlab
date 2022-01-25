function twoP_plotAUC
%%

S = twoP_settings;
colAnimalID = 1; colDate = 2; colLocation = 3; colDepth = 4; colExpertise = 5; colSession = 6;
aucFileName = 'AUC_6.mat';
filelist = dir(fullfile(S.dir.imagingRootDir,['**\' aucFileName]));
folderList = {filelist.folder};

splitCells = cellfun(@(x) regexp(x,filesep,'split'),folderList,'UniformOutput',false);
animalID = cellfun(@(x) x{7},splitCells,'UniformOutput',false); % Parsed animalID from the directory list
sessionID = cellfun(@(x) x{9},splitCells,'UniformOutput',false); % parsed sessionID from the directory list
idxExps = cell2mat(cellfun(@(x,y) find(strcmp(S.exps(:,1),x) & strcmp(S.exps(:,6),y)), animalID, sessionID, 'UniformOutput',false));
subExps = [S.exps(idxExps,1:6) S.depth(idxExps)];
cAnimal = subExps(:,1); cDate = subExps(:,2); cLocation = subExps(:,3); cDepthNumber = subExps(:,4);
cExpertise = subExps(:,5); cSession = subExps(:,6); cDepth = subExps(:,7);

%% Load and combine AUC of select sessions
close all

sAnimal = {'CSP27';'CSP30'};
sExpertise = {'Expert'};
sDepth = {'Intermediate'};
sLocation = {'ALM'};

sIdx = find(ismember(cAnimal,sAnimal) & ...
    strcmp(cExpertise,sExpertise) & ...
    strcmp(cDepth,sDepth) & ...
    strcmp(cLocation,sLocation));
sMeta = subExps(sIdx,:);
allAUC = cell(length(sIdx),1);
shufAUC = cell(length(sIdx),1);
idxRed = cell(length(sIdx),1);

parfor i = 1:length(sIdx)
    ii=sIdx(i);
    AUC = load(fullfile(filelist(ii).folder,filelist(ii).name));
    AUC = AUC.AUC;
    idxCell = readNPY(fullfile(filelist(ii).folder,'iscell.npy'));
    idxRedTemp = readNPY(fullfile(filelist(ii).folder,'redcell.npy'));
    idxRedTemp = logical(idxRedTemp(logical(idxCell(:,1))));
    allAUC{i} = AUC.allAUC;
    shufAUC{i} = AUC.shufAUC;
    idxRed{i} = idxRedTemp;
end

allAUC = cell2mat(allAUC); 
shufAUC = cell2mat(shufAUC);
idxRed = cell2mat(idxRed);
idxCellTypes = twoP_indexMatrix(idxRed);

%% Find AUCs that are significantly different from the mean

sigAUC = nan(size(allAUC,1),size(allAUC,2));
uShuf = mean(shufAUC,3)+2*std(shufAUC,0,3); lShuf = mean(shufAUC,3)-2*std(shufAUC,0,3);
uAUC = (allAUC >= uShuf).*allAUC; lAUC = (allAUC <= lShuf).*allAUC;
sigAUC = uAUC + lAUC; sigAUC(sigAUC==0) = nan;

for j = 1:size(sigAUC,2)
    for jj = 1:size(idxCellTypes,2)
        side{1}(j,jj) = sum(sigAUC(idxCellTypes(:,jj),j) > 0.5);
        side{2}(j,jj) = sum(sigAUC(idxCellTypes(:,jj),j) < 0.5);
    end
end

%% Plot AUC fractions of cell types in bar graphs

hFigBars = figure('Position',[500, 500, 600, 300]);
hTilesBars = tiledlayout(1,numel(side));
bar_colors = {'black','red','green'};
x_tick_text = {'Init','EarlyStim','LateStim','Delay','EarlyResp','LateResp'};
legend_text = {'\color{black} All','\color{red} tdT+','\color{green} tdT-'};
tile_title = {'Contra','Ipsi'};

for i = 1:length(side)
        ax{i} = nexttile;
        hold on
        colororder(bar_colors)
        b = bar(side{i}./sum(idxCellTypes),'EdgeColor','none');
        title(tile_title{i});
        y_lim_max(i) = max(get(gca,'ylim'));
        xticklabels(x_tick_text);
        hL = legend(legend_text,...
            'Location','northwest',...
            'Box','off');
        if i == 1
            ylabel('Fraction of cells (%)');
        end
end
linkaxes([ax{:}],'y');
ylim([0 max(y_lim_max)]);
fig_configAxis(gca);
title(hTilesBars,[cell2mat(sAnimal') ' ' sExpertise{:} ' ' sLocation{:} ' ' sDepth{:} ' depth']);
exportgraphics(hFigBars,fullfile(S.dir.imagingRootDir,'AUC',[cell2mat(sAnimal') '_' sExpertise{:} '_' sLocation{:} '_' sDepth{:} '_AUC_fractions_bars.pdf']));


%% Plot histograms of combined/concatenated AUC data

[hFigure, hTiles] = twoP_histEpoches(sigAUC,idxCellTypes,20);
title(hTiles,[cell2mat(sAnimal') ' ' sExpertise{:} ' ' sLocation{:} ' ' sDepth{:} ' depth']);
exportgraphics(hFigure,fullfile(S.dir.imagingRootDir,'AUC',[cell2mat(sAnimal') '_' sExpertise{:} '_' sLocation{:} '_' sDepth{:} '_AUC_combined.pdf']));

    
%% Load and plot selected animals

animalSel = {'CSP','Fez','Plex'};
idxAnimal = find(contains({filelist.folder},animalSel));

for i = 1:length(idxAnimal)
    %%
    close all
    
    animal = cAnimal{i};
    session = cSession{i};
    location = cLocation{i};
    depth = cDepthNumber{i};
    
    load(fullfile(filelist(i).folder,filelist(i).name));
    idxCell = readNPY(fullfile(filelist(i).folder,'iscell.npy'));
    idxRed = readNPY(fullfile(filelist(i).folder,'redcell.npy'));
    idxRed = logical(idxRed(logical(idxCell(:,1))));

    numBins = 30;
    limX = [0 1];
    sigAUC = nan(size(AUC.allAUC,1),size(AUC.allAUC,2));
    binCounts = zeros(size(AUC.allAUC,2),numBins);
    uShuf = mean(AUC.shufAUC,3)+2*std(AUC.shufAUC,0,3);
    lShuf = mean(AUC.shufAUC,3)-2*std(AUC.shufAUC,0,3);
    uAUC = (AUC.allAUC >= uShuf).*AUC.allAUC;
    lAUC = (AUC.allAUC <= lShuf).*AUC.allAUC;
    sigAUC = uAUC + lAUC; sigAUC(sigAUC==0) = nan;

    hFigure = figure('Position',[500 500 1200 250]);
    hTiles = tiledlayout(1,size(AUC.allAUC,2));
    for j = 1:size(sigAUC,2)
        binCounts(j,:) = histcounts(sigAUC(:,j), 30,'BinLimits',limX);
    end
%     limY = [min(binCounts(:)) ceil(max(binCounts(:))/10)*10];
    for j = 1:size(sigAUC,2)
        nAll(i,j) = round(sum(sigAUC(:,j)>0.5)/sum(sigAUC(:,j)<0.5),2,'significant');
        nR(i,j) = round(sum(sigAUC(idxRed,j)>0.5)/sum(sigAUC(idxRed,j)<0.5),2,'significant');
        nNR(i,j) = round(sum(sigAUC(~idxRed,j)>0.5)/sum(sigAUC(~idxRed,j)<0.5),2,'significant');

        nexttile
        hold on
        histogram(sigAUC(~isnan(sigAUC(:,j)),j), numBins,...
            'BinLimits',limX,...
            'DisplayStyle','stairs',...
            'EdgeColor','k',...
            'LineWidth',1,...
            'Normalization','probability');
        histogram(sigAUC(~isnan(sigAUC(idxRed,j)),j), numBins,...
            'BinLimits',limX,...
            'DisplayStyle','stairs',...
            'EdgeColor','r',...
            'LineWidth',1,...
            'Normalization','probability');
        histogram(sigAUC(~isnan(sigAUC(~idxRed,j)),j), numBins,...
            'BinLimits',limX,...
            'DisplayStyle','stairs',...
            'EdgeColor','g',...
            'LineWidth',1,...
            'Normalization','probability');
        hL = legend(['All: ' num2str(sum(sigAUC(:,j)>0.5)) '/' num2str(sum(sigAUC(:,j)<0.5)) '=' num2str(nAll(i,j))],...
            ['tdT+: '  num2str(sum(sigAUC(idxRed,j)>0.5)) '/' num2str(sum(sigAUC(idxRed,j)<0.5)) '=' num2str(nR(i,j))],...
            ['tdT-: '  num2str(sum(sigAUC(~idxRed,j)>0.5)) '/' num2str(sum(sigAUC(~idxRed,j)<0.5)) '=' num2str(nNR(i,j))],...
            'Location','northwest',...
            'Box','off');
        offsetAxes(gca);
        fig_configAxis(gca);
    end
    title(hTiles,[animal ' ' session ' ' location ' ' depth '\mum']);
    exportgraphics(hFigure,fullfile(S.dir.imagingRootDir,'AUC',[animal '_' session '_AUC.pdf']));
end

%% Plot the ratio of ipsi/contra AUC neuron counts

close all;

plotAnimal = {'CSP30'};
plotLocation = 'MM';
plotDepth = 'Intermediate';

markerstyles = {'o', '^', '*'}; 
markercolors = {'r', 'g', 'k'};
markerlegends = {'tdT+', 'tdT-', 'All'};

idxSummary = ismember(subExps(:,1),plotAnimal) & ...
    strcmp(subExps(:,3),plotLocation) & ...
    strcmp(subExps(:,7),plotDepth);
xCell = {datenum(subExps(idxSummary,2))};
yCell = {nR(idxSummary,4),nNR(idxSummary,4),nAll(idxSummary,4)};

if length(xCell) == 1
    xCell = repmat(xCell,length(yCell),1);
elseif length(xCell) ~= length(yCell)
    disp('The number of elements between the x-axis and the y-axis are different! Please check again.');
    return
end

hF = figure; hold on;
set(hF,'Position',[500 500 500 250]);

for i = 1:length(yCell)
hP = plot(xCell{i},yCell{i}); 

set(hP, 'color',markercolors{i},...
    'linestyle','none',...
    'marker', markerstyles{i},...
    'MarkerSize',5);

pHandles(i) = hP;
end

xlabel('Session date');
ylabel('Ratio: ipsi to contra neurons');
title([cell2mat(plotAnimal) ' ' plotLocation ' ' plotDepth]);

hL = legend(pHandles,markerlegends,'Box','off');
lpos = get(hL, 'position');
lpos(1) = lpos(1) + 0.10; 
lpos(2) = lpos(2) - 0.02;
set(hL, 'location', 'northeast', 'box', 'on');

axis tight; ylims = get(gca, 'ylim');
if max(ylims) >= 5
    ylim([0 max(abs(ylims))/2]);
else
    ylim([0 max(abs(ylims))]);
end

xtickangle(gca,90);
datetick('x','mmm-dd');
fig_configAxis(gca); 
offsetAxes(gca);
exportgraphics(hF,fullfile(S.dir.imagingRootDir,'AUC',[cell2mat(plotAnimal) '_' plotLocation '_' plotDepth '_AUC_ratio_neuron_counts.pdf']));
