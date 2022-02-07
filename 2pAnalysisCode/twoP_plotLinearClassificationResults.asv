function twoP_plotAUC
%%

S = twoP_settings;
colAnimalID = 1; colDate = 2; colLocation = 3; colDepth = 4; colExpertise = 5; colSession = 6;
aucFileName = 'LR_136.mat';

filelist = dir(fullfile(S.dir.imagingRootDir,['**\' aucFileName]));
folderList = {filelist.folder};

splitCells = cellfun(@(x) regexp(x,filesep,'split'),folderList,'UniformOutput',false);
animalID = cellfun(@(x) x{7},splitCells,'UniformOutput',false); % Parsed animalID from the directory list
sessionID = cellfun(@(x) x{9},splitCells,'UniformOutput',false); % parsed sessionID from the directory list
idxExps = cell2mat(cellfun(@(x,y) find(strcmp(S.exps(:,1),x) & strcmp(S.exps(:,6),y)), animalID, sessionID, 'UniformOutput',false));
subExps = [S.exps(idxExps,1:6) S.depth(idxExps)];
cAnimal = subExps(:,1); cDate = subExps(:,2); cLocation = subExps(:,3); cDepthNumber = subExps(:,4);
cExpertise = subExps(:,5); cSession = subExps(:,6); cDepth = subExps(:,7);

%% Loads logistic regression analysis, takes the mean of an epoch, and combines multiple sessions into rows of a cell

close all

sAnimal = {''};
sExpertise = {''};
sDepth = {''};
sLocation = {''};

if ismissing(sAnimal)
    sAnimal = unique(cAnimal);
end

if ismissing(sExpertise)
    sExpertise = unique(cExpertise);
end

if ismissing(sDepth)
    sDepth = unique(cDepth);
end

if ismissing(sLocation)
    sLocation = unique(cLocation);
end


sIdx = find(contains(cAnimal,sAnimal) & ...
    contains(cExpertise,sExpertise) & ...
    contains(cDepth,sDepth) & ...
    contains(cLocation,sLocation));

% This defines which epoch is averaged
eIdx = [S.segFrames(3) S.segFrames(4)-1];

A = cell(length(sIdx),1);
betaR = cell(length(sIdx),1);
betaU = cell(length(sIdx),1);
shufA = cell(length(sIdx),1);
R = cell(length(sIdx),1);
U = cell(length(sIdx),1);
UR = cell(length(sIdx),1);
UU = cell(length(sIdx),1);

sMeta = subExps(sIdx,:);

folderList = {filelist.folder};
filenameList = {filelist.name};
parfor i = 1:length(sIdx)
    try
    ii=sIdx(i);
    lr = load(fullfile(folderList{ii},filenameList{ii}),'lr'); lr = lr.lr;
    twoP_plotSingleSessionLinearClassification(sMeta(i,:),lr)
    idxCell = readNPY(fullfile(filelist(ii).folder,'iscell.npy'));
    idxRedTemp = readNPY(fullfile(filelist(ii).folder,'redcell.npy'));
    idxRed{i} = logical(idxRedTemp(logical(idxCell(:,1))));
    A{i} = mean(lr.cvAcc(eIdx(1):eIdx(2)),2,'omitnan');
    betaR{i} = mean(lr.bMaps(idxRed{i},eIdx(1):eIdx(2)),2,'omitnan');
    betaU{i} = mean(lr.bMaps(~idxRed{i},eIdx(1):eIdx(2)),2,'omitnan');
    shufA{i} = mean(lr.cvAccShuf(eIdx(1):eIdx(2)),2,'omitdbquitnan');
    R{i} = mean(lr.cvAccRed(eIdx(1):eIdx(2)),2,'omitnan');
    U{i} = mean(lr.cvAccNR(:,eIdx(1):eIdx(2)),2,'omitnan');
    UR{i} = mean(lr.cvAccMixedUR(:,eIdx(1):eIdx(2)),2,'omitnan');
    UU{i} = mean(lr.cvAccMixedUU(:,eIdx(1):eIdx(2)),2,'omitnan');
    disp(['Loaded ' fullfile(filelist(ii).folder,filelist(ii).name)]);
    end
end

idxRedAll = vertcat(idxRed{:});
idxCellTypes = twoP_indexMatrix(vertcat(idxRed{:}));


%%
time_vec = 1:length(lr.cvAcc);
hFig = figure; hold on;
boundedline(time_vec,mean(lr.cvAccNR),std(lr.cvAccNR,0,1,'omitnan')./sqrt(size(lr.cvAccNR,1)),'g','nan','gap',...
    'transparency',0.1);
boundedline(time_vec,mean(lr.cvAccMixedUR),std(lr.cvAccMixedUR,0,1,'omitnan')./sqrt(size(lr.cvAccMixedUR,1)),'m','nan','gap',...
    'transparency',0.1);
boundedline(time_vec,mean(lr.cvAccMixedUU),std(lr.cvAccMixedUU,0,1,'omitnan')./sqrt(size(lr.cvAccMixedUU,1)),'b','nan','gap',...
    'transparency',0.1)
line(time_vec,lr.cvAcc,'color','k');
line(time_vec,lr.cvAccShuf,'color',[0.5 0.5 0.5]);
line(time_vec,lr.cvAccRed,'color','r');
% str_title_label = strjoin(string(horzcat(sMeta(i,1),sMeta(i,6),sMeta(i,3),sMeta(i,5),sMeta(i,end))));
title(str_title_label)
% ytickformat(gca, '%g%%');
ax = gca;
set(ax,'ytick',0.3:0.1:1);
offsetAxes(gca);
fig_configAxis(gca);
% exportgraphics(hF,fullfile(S.dir.imagingRootDir,'LogisticRegression',['AllSessionsLogisticRegression_' sAnimal{:} '.pdf']));

%%
xval = 1:size(sIdx,1);
x_labels = horzcat(sMeta(:,1),sMeta(:,6),sMeta(:,3),sMeta(:,5),sMeta(:,end));
str_x_labels = string(x_labels);
for i = 1:size(str_x_labels,1)
    x_tick_labels(i) = strjoin(str_x_labels(i,:));
end

hF = figure('position',[500 500 1500 500]); hold on;
h = ploterr(xval,cell2mat(A),[],[],'sk');
hShuf = ploterr(xval,cell2mat(shufA),[],[]);
set(hShuf,'marker','s',...
    'linestyle','none',...
    'markeredgecolor',[0.5 0.5 0.5]);
hR = ploterr(xval,cell2mat(R),[],[],'*r');
hU = ploterr(xval,cellfun(@mean,U),[],cellfun(@std,U),'og');
hUR = ploterr(xval,cellfun(@mean,UR),[],cellfun(@std,UR),'*m');
hUU = ploterr(xval,cellfun(@mean,UU),[],cellfun(@std,UU),'ob');
xticks(xval);
xticklabels(x_tick_labels)
l = legend('All Cells','All Cells - Shuffled',...
    'tdT+','','tdT- subsampled', '',...
    'Mixed','','Mixed subsampled');
set(l,'location','northwest',...
    'box','off');
ylabel('Classifier accuracy (%)');
title(['Classifier accuracy during delay epoch across sessions - ' sAnimal{:}])
% set(hUR,'marker',)
offsetAxes(gca);
fig_configAxis(gca);
exportgraphics(hF,fullfile(S.dir.imagingRootDir,'LogisticRegression',['AllSessionsLogisticRegression_' sAnimal{:} '.pdf']));

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
