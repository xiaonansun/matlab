function twoP_plotLinearClassificationResults
%% 2022-02-21: at this time, please run one cell at a time. 

S = twoP_settings;
colAnimalID = 1; colDate = 2; colLocation = 3; colDepth = 4; colExpertise = 5; colSession = 6;
findFileName = 'LR_6.mat';

filelist = dir(fullfile(S.dir.imagingRootDir,['**\' findFileName]));
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

% This defines the averaged epoch
numTimePoints = str2num(cell2mat(regexp(findFileName,'\d*','Match')));
if numTimePoints > 6
eIdx = [S.segFrames(3) S.segFrames(4)-1];
else
    eIdx = [4 4];
end

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
%     twoP_plotSingleSessionLinearClassification(lr,sMeta(i,:));
    idxCell = readNPY(fullfile(filelist(ii).folder,'iscell.npy'));
    idxRedTemp = readNPY(fullfile(filelist(ii).folder,'redcell.npy'));
    idxRed{i} = logical(idxRedTemp(logical(idxCell(:,1))));
    A{i} = mean(lr.cvAcc(eIdx(1):eIdx(2)),2,'omitnan');
    betaR{i} = mean(lr.bMaps(idxRed{i},eIdx(1):eIdx(2)),2,'omitnan');
    betaU{i} = mean(lr.bMaps(~idxRed{i},eIdx(1):eIdx(2)),2,'omitnan');
    shufA{i} = mean(lr.cvAccShuf(:,eIdx(1):eIdx(2)),2,'omitnan'); 
    Red{i} = mean(lr.cvAccRed(eIdx(1):eIdx(2)),2,'omitnan');
    R{i} = mean(lr.cvAccR(:,eIdx(1):eIdx(2)),2,'omitnan');
    U{i} = mean(lr.cvAccU(:,eIdx(1):eIdx(2)),2,'omitnan');
    UR{i} = mean(lr.cvAccMixedUR(:,eIdx(1):eIdx(2)),2,'omitnan');
    UU{i} = mean(lr.cvAccMixedUU(:,eIdx(1):eIdx(2)),2,'omitnan');
    disp(['Loaded ' fullfile(filelist(ii).folder,filelist(ii).name)]);
    end
end

idxRedAll = vertcat(idxRed{:});
idxCellTypes = twoP_indexMatrix(vertcat(idxRed{:}));

%% Plot histogram of beta values
close all

sel = {'CSP','Intermediate','Expert','ALM'};
% sel = {'CSP','Plex','Fez','Deep','Intermediate','Superficial','Trained','Expert','ALM','MM'};
threshBeta = 0.01;
subMeta = sMeta(:,[1,3,5,7]);
subIdx = find(prod(contains(subMeta,sel),2));

allBetaU = vertcat(betaU{subIdx});
allBetaU(abs(allBetaU)>threshBeta)=nan;
allBetaR = vertcat(betaR{subIdx});
allBetaR(abs(allBetaR)>threshBeta)=nan;
catU = repmat({'U'},size(allBetaU,1),size(allBetaU,2));
catR = repmat({'R'},size(allBetaR,1),size(allBetaR,2));
catAll = vertcat(catU,catR);
allBeta = vertcat(allBetaU,allBetaR);

x_labels = horzcat(sMeta(:,1),sMeta(:,6));
str_x_labels = string(x_labels);
for i = 1:size(str_x_labels,1)
    x_tick_labels(i) = strjoin(str_x_labels(i,:));
end
color_labeled = {'r','m'}; color_unlabeled = {'g','b'}; 

xval = 1:length(subIdx);
cA = A; cRed = Red; cR = R; cU = U; cUU = UU; cUR = UR; cshufA = shufA;
cA(logical(cellfun(@isempty,cA))) = {single(nan)}; cA = cell2mat(cA);
cRed(logical(cellfun(@isempty,cRed))) = {single(nan)}; cRed = cell2mat(cRed);
cshufA(logical(cellfun(@isempty,cshufA))) = {single(nan)}; 
cR(logical(cellfun(@isempty,cR))) = {single(nan)};
cU(logical(cellfun(@isempty,cU))) = {single(nan)};
cUU(logical(cellfun(@isempty,cUU))) = {single(nan)};
cUR(logical(cellfun(@isempty,cUR))) = {single(nan)};

% Hypthesis testing
[~,p_Welch] = ttest2(allBetaU,allBetaR); %Welch's t-test
p_KW=kruskalwallis(allBeta,catAll,'off'); % Kruskal-Wallis test

% if ishandle(hFigHist)
%     figure(hFigHist);
% else
    hFigHist = figure('position',[500 500 350 300]);
% end

hold on;
hHist(1) = histogram(allBetaU,100,'FaceColor','g'); 
hHist(2) = histogram(allBetaR,100,'FaceColor','r'); 
for i = 1:length(hHist)
hHist(i).EdgeColor = 'none';
end
ax = gca; xlabel('\beta-weight'); ylabel('Number of neurons');
% ax.Units = 'pixels';
hLeg = legend(['\color{green}tdT- (n='  num2str(length(allBetaU)) ')'],['\color{red}tdT+ (n='  num2str(length(allBetaR)) ')'],...
    'box','off',...
    'color','none');

offsetAxes(gca);
fig_configAxis(gca);
title(strjoin(string(sel)));
text(ax,ax.XTick(end),ax.YTick(4),['\itp \rm= ' num2str(round(p_Welch,3,'significant'))],...
    'HorizontalAlignment','right');
text(ax,ax.XTick(end),mean(ax.YTick(3:4)),['\itp \rm= ' num2str(round(p_KW,3,'significant'))],...
    'HorizontalAlignment','right');

insetAx = axes('Units','normalized',...
    'Position',[ax.Position(1)*1.2 ax.Position(2)*4 ax.Position(3)*0.25 ax.Position(4)*0.25]); hold on;
histogram(insetAx, allBetaU,100,'FaceColor','none',...
    'EdgeColor','g',...
    'LineWidth',1,...
    'DisplayStyle','stairs',...
    'Normalization','probability');
histogram(insetAx, allBetaR,100,'FaceColor','none',...
    'EdgeColor','r',...
    'LineWidth',1,...
    'DisplayStyle','stairs',...
    'Normalization','probability');
insetAx.Box = 'off'; insetAx.Color = 'none';
saveHistFigurePath = fullfile(S.dir.imagingRootDir, 'LogisticRegression', [horzcat(sel{:}) '_histogram.pdf']);
exportgraphics(hFigHist,saveHistFigurePath); disp(['Figure saved as ' saveHistFigurePath]);

hFigAcc = figure('Position',[500 500 1000 300]);
tiledlayout(1,5)
tPure = nexttile(1,[1 2]); hold on;
hA(1) = ploterr(xval,cA(subIdx),[],[],'sk'); leg1(1) = hA(1); 
hShuf = ploterr(xval+0.1,cellfun(@mean,cshufA(subIdx)),[],cellfun(@std,cshufA(subIdx))./sqrt(cellfun(@length,cshufA(subIdx))),'sk',...
    'abshhxy',0); leg1(2) = hShuf(1); marShuf(1) = hShuf(1);
hR = ploterr(xval+0.2,cellfun(@mean,R(subIdx)),[],cellfun(@std,R(subIdx))./sqrt(cellfun(@length,R(subIdx))),'^r',...
    'abshhxy',0); leg1(3) = hR(1); marLabeled(1) = hR(1);
hU = ploterr(xval+0.3,cellfun(@mean,U(subIdx)),[],cellfun(@std,U(subIdx))./sqrt(cellfun(@length,U(subIdx))),'og',...
    'abshhxy',0); leg1(4) = hU(1); marUnlabeled(1) = hU(1);
ax(1) = gca;
l(1) = legend(leg1(:), 'All Cells','All Cells - Shuf',...
    'tdT+','tdT- subsamp');

tMixed = nexttile(3,[1 2]); hold on;
hA(2) = ploterr(xval,cA(subIdx),[],[],'sk'); leg2(1) = hA(2);
hShuf = ploterr(xval+0.1,cellfun(@mean,cshufA(subIdx)),[],cellfun(@std,cshufA(subIdx))./sqrt(cellfun(@length,cshufA(subIdx))),'sk',...
    'abshhxy',0); leg2(2) = hShuf(1); marShuf(2) = hShuf(1);
hUR = ploterr(xval+0.2,cellfun(@mean,UR(subIdx)),[],cellfun(@std,UR(subIdx))./sqrt(cellfun(@length,UR(subIdx))),'^m',...
    'abshhxy',0); leg2(3) = hUR(1); marLabeled(2) = hUR(1);
hUU = ploterr(xval+0.3,cellfun(@mean,UU(subIdx)),[],cellfun(@std,UU(subIdx))./sqrt(cellfun(@length,UU(subIdx))),'ob',...
    'abshhxy',0); leg2(4) = hUU(1); marUnlabeled(2) = hUU(1);
ax(2) = gca;
l(2) = legend(leg2(:),'All Cells','All Cells - Shuf', ...
    'Mixed','Mixed subsampled');

for i = 1:length(ax)
    offsetAxes(ax(i));
    fig_configAxis(ax(i));
    hA(i).MarkerFaceColor = 'k'; hA(i).MarkerEdgeColor = 'w';
    marShuf(i).MarkerFaceColor = [0.5 0.5 0.5]; marShuf(i).MarkerEdgeColor = 'w';
    marLabeled(i).MarkerFaceColor = color_labeled{i};
    marLabeled(i).MarkerEdgeColor = 'w';
    marUnlabeled(i).MarkerFaceColor = color_unlabeled{i};
    marUnlabeled(i).MarkerEdgeColor = 'w';
    for j = 1:length(ax(i).Children)
        if ax(i).Children(j).LineStyle == '-'
            ax(i).Children(j).LineWidth = 1;
        end
    end
    ax(i).XTick = xval;
    ax(i).XTickLabel = cellstr(x_tick_labels(subIdx));
    ax(i).YLim = [0.4 1];
    ax(i).YLabel.String = 'Classifier accuracy (%)';
    ax(i).Title.String = 'Session breakdown';
    l(i).Location = 'northwest';
    l(i).Box = 'off';
end

tSummary = nexttile(5); hold on; % don't need to use ax = gca if using tile handle. Tile handle replaces axes handle.
legend_labels = {'All','All (Shuf)','tdT+','tdT-','Mixed','tdT- (2X)'};

marker_color = {'k',[0.5 0.5 0.5],'r','g','m','b'};
marker_color = flip(marker_color);
hA = ploterr(1,mean(cA(subIdx),'omitnan'),[],std(cA(subIdx),0,'omitnan'),'sk',...
    'abshhxy',0);
hShuf = ploterr(1+0.1,mean(cellfun(@mean,cshufA(subIdx)),'omitnan'),[],std(cellfun(@mean,cshufA(subIdx)),0,'omitnan')./sqrt(length(cellfun(@mean,cshufA(subIdx)))),'sk',...
    'abshhxy',0);
hR = ploterr(1+0.2,mean(cellfun(@mean,R(subIdx)),'omitnan'),[],std(cellfun(@mean,R(subIdx)),0,'omitnan')./sqrt(length(cellfun(@mean,R(subIdx)))),'^r',...
    'abshhxy',0);
hU = ploterr(1+0.3,mean(cellfun(@mean,U(subIdx)),'omitnan'),[],std(cellfun(@mean,U(subIdx)),0,'omitnan')./sqrt(length(cellfun(@mean,U(subIdx)))),'og',...
    'abshhxy',0);
hUR = ploterr(1+0.4,mean(cellfun(@mean,UR(subIdx)),'omitnan'),[],std(cellfun(@mean,UR(subIdx)),0,'omitnan')./sqrt(length(cellfun(@mean,UR(subIdx)))),'^m',...
    'abshhxy',0); 
hUU = ploterr(1+0.5,mean(cellfun(@mean,UU(subIdx)),'omitnan'),[],std(cellfun(@mean,UU(subIdx)),0,'omitnan')./sqrt(length(cellfun(@mean,UU(subIdx)))),'ob',...
    'abshhxy',0); 
tSummary.XAxis.Visible = 'off';
cntLine = 1; cntMarker = 1;
for i = 1:length(tSummary.Children)
    if tSummary.Children(i).LineStyle == '-'
        tSummary.Children(i).Color = marker_color{cntLine};
%         hText(cntLine) = text(min(tSummary.Children(i).XData), ...
%             min(tSummary.Children(i).YData), ...
%             legend_labels{cntLine},...
%             'Rotation',90,...
%             'VerticalAlignment','top',...
%             'Color',marker_color{cntLine});
        cntLine = cntLine +1;
    end
    if ~strcmp(tSummary.Children(i).Marker,'none')
        tSummary.Children(i).MarkerFaceColor = marker_color{cntMarker};
        tSummary.Children(i).MarkerEdgeColor = 'w';
        cntMarker = cntMarker + 1;
    end
end
offsetAxes(tSummary);
fig_configAxis(tSummary);
tSummary.YLim = [0.4 1];
[~,fn,~] = fileparts(findFileName); 
saveAccFigurePath = fullfile(S.dir.imagingRootDir, 'LogisticRegression', [horzcat(sel{:}) '_' fn '_cvAcc.pdf']);
exportgraphics(hFigAcc,saveAccFigurePath); disp(['Figure saved as ' saveAccFigurePath]);

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
