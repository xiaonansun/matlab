

%%
% toLoad = {'CSP','Trained','MM','Superficial';
%     'CSP','Trained','MM','Intermediate';
%     'CSP','Trained','MM','Deep'};

toLoad = {'CSP','Trained','ALM','Superficial';
    'CSP','Trained','ALM','Intermediate';
    'CSP','Trained','ALM','Deep'};

% toLoad = {'Fez','Trained','ALM','Deep';
%     'Fez','Trained','MM','Deep'};

for i = 1:size(toLoad,1)
    [E{i},M{i}]= twoP_loadSessions(toLoad{i,1},toLoad{i,2},toLoad{i,3},toLoad{i,4});
end
% E: all data; M: meta-data

fName = [];
for i = 1:size(toLoad,1)
    fName = [fName toLoad{i,1} toLoad{i,3} toLoad{i,4}];
    if i < size(toLoad,1)
        fName = [fName '_vs_'];
    end
end

%%
S = twoP_settings;
sPerFrame = S.msPerFrame/1000;
trialType = 4; %1: stim, 2: response; 3: rewarded; 4: error

fData = length(E);
idx = find(~cellfun(@isempty, E{1}(:,4)));
trialTypeNames = fieldnames(E{1}{idx(1),5}.sub.names);
fldBhv = fieldnames(E{1}{idx(1),5}.(trialTypeNames{trialType}));
posNeurons = []; negneurons = []; diffNeurons = [];

for f = 1:length(E)
    allData = E{f};
    idxRow = find(~cellfun(@isempty, allData(:,4)));
    for j = 1:length(fldBhv)
        nVc = [];
        idxRed = [];
        countNan = [];
        for i = 1:length(idxRow)
            tempVc = mean(allData{idxRow(i),4}(:,:,allData{idxRow(i),5}.(trialTypeNames{trialType}).(fldBhv{j})),3,'omitnan');
            tempVc(tempVc > 3) = NaN;
            nVc = [nVc ; tempVc];
            idxRed = [idxRed;allData{idxRow(i),6}];
            countNan = [countNan;allData{idxRow(i),7}];
        end
        cellCount(f).notRed = sum(~idxRed); cellCount(f).Red = sum(idxRed);
        allNan = sum(countNan);
        filterNan = ones(1,length(allNan));
        [valPks,idxPks] = findpeaks(allNan);
        mu = 0.5;
        for i = 1:length(idxPks)
            [valMin(i),idxMin(i)]=min(abs(allNan(S.segFrames(i):idxPks(i))-mu*valPks(i)));
            filterNan(S.segFrames(i)+idxMin(i):idxPks(i))=nan;
        end
        t = 1:size(nVc,2);
        posNeurons(j,:,f) = mean(nVc(find(idxRed),:),1,'omitnan') .* filterNan;
        negNeurons(j,:,f) = mean(nVc(find(~idxRed),:),1,'omitnan') .* filterNan;
        diffNeurons(j,:,f) = posNeurons(j,:,f)-negNeurons(j,:,f);
    end
end
%% Plot figures
close all;

epochAligned = nan(size(diffNeurons,1),size(diffNeurons,2),size(diffNeurons,3));
segFrames = [S.segFrames(1:4)+1 S.segFrames(5)];
for j = 1:size(diffNeurons,3)
    for i = 1:length(segFrames)-1
        epochAligned(:,segFrames(i):segFrames(i+1)-1,j) = ...
            diffNeurons(:,segFrames(i):segFrames(i+1)-1,j) - diffNeurons(:,segFrames(i),j);
    end
end

figEpoch = figure(2);
set(figEpoch,'Position',[100 100 400 300*size(figEpoch,3)]);
tEpoch = tiledlayout(length(E),length(fldBhv),...
    'TileSpacing','tight',...
    'Padding','tight');
YLim = squeeze([min(epochAligned,[],[1 2],'omitnan') max(epochAligned,[],[1 2],'omitnan')])';
t = 0:sPerFrame:sPerFrame*(size(epochAligned,2)-1);

hold on;
clrOrange = [0.9290 0.6940 0.1250]; clrMagenta = [1 0 1];
if size(epochAligned,1) == 3
    lineColor = [[0 0 0]; clrMagenta; clrOrange];
    lineMarker = {'-';'-';'-'};
elseif size(epochAligned,1) == 4
    lineColor = [clrMagenta; clrOrange; clrMagenta; clrOrange];
    lineMarker = {':';':';'-';'-'};
end

for i = 1:size(epochAligned,3)
    nexttile
    for j = 1:size(epochAligned,1)
        line(t,epochAligned(j,:,i),...
            'LineStyle',lineMarker{j},...
            'Color',lineColor(j,:),...
            'LineWidth',2);
        ax = gca;
        ax = fig_configAxis(ax);
        if i ~= size(epochAligned,3); ax.XTickLabel = []; ax.XAxis.Color = 'none'; end
    end
    xlim([min(t) max(t)]);
    title({[M{i}.cellType ' ' M{i}.location ' ' M{i}.depth ' ' trialTypeNames{trialType} ' trials']});
    ylabel({'Z_{diff}'});
    
    if i == 2
        legend(allData{idxRow(i),5}.sub.names.(trialTypeNames{trialType}),...
            'Location','Northwest',...
            'Box','off');
    end
    if i == 3; xlabel('Time (s)'); end
end

exportgraphics(figEpoch,fullfile(S.dir.imagingRootDir,[fName '_' trialTypeNames{trialType} '_Zdiff.pdf']));


% Plot and save figure(s)

figAllPlots = figure(1);
set(figAllPlots,'Position',[100 100 1500 800]);
tFig = tiledlayout(length(E),length(fldBhv),...
    'TileSpacing','compact',...
    'Padding','compact');
leftYLim = squeeze([min([posNeurons;negNeurons],[],[1 2],'omitnan') max([posNeurons;negNeurons],[],[1 2],'omitnan')])';
rightYLim = squeeze([min(diffNeurons,[],[1 2],'omitnan') max(diffNeurons,[],[1 2],'omitnan')])';
t = 0:sPerFrame:sPerFrame*(size(posNeurons,2)-1);

for f = 1:length(E)
    for j = 1:size(posNeurons,1)
        nexttile
        
        hold on;
        
        yyaxis left
        line(t,posNeurons(j,:,f),...
            'Color','r',...
            'LineWidth',1);
        line(t,negNeurons(j,:,f),...
            'Color','g',...
            'LineWidth',1);
        ax1 = gca;
        if j > 1; ax1.YTickLabel = []; ax1.YAxis(1).Visible = 'off'; end
        if j == 1; ylabel({[M{f}.cellType ' ' M{f}.location]; M{f}.depth}); end
        ylim(leftYLim(f,:));
        
        yyaxis right
        line(t,diffNeurons(j,:,f),...
            'Color','k',...
            'LineWidth',3);
        ylim(rightYLim(f,:));
        xline(S.segFrames(1));
        ax2 = gca;
        ax2 = fig_configAxis(ax2);
        
        if f == 1; title(allData{idxRow(i),5}.sub.names.(trialTypeNames{trialType}){j}); end
        if f == length(E); xlabel('Time (s)'); end
        if j == 1
            legend(['tdT_+ (n=' num2str(cellCount(f).Red) ')'],...
                ['tdT_- (n=' num2str(cellCount(f).notRed) ')'],...
                'tdT_+ - tdT_-',...
                'Location','Northwest',...
                'Box','off'); 
        end
        if f ~= length(E); ax2.XTickLabel = []; ax2.XAxis.Color = 'none'; end
        if j == 4; ylabel('Z-score_+ - Z-score_-'); end
        if j < 4
            ax2.YTickLabel = []; ax2.YAxis(2).Visible = 'off';
        else
            ax2.YAxis(2).Color = [0 0 0];
        end
        
        
        xlim([min(t) max(t)])
    end
end

exportgraphics(figAllPlots,fullfile(S.dir.imagingRootDir,[fName '_' trialTypeNames{trialType} '.pdf']));