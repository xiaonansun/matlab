

%% Load data files

clear E M
% toLoad = {'CSP','Trained','ALM','Superficial';...
%     'CSP','Trained','ALM','Intermediate';...
%     'CSP','Trained','ALM','Deep'};

% toLoad = {'CSP','Trained','MM','Superficial';
%     'CSP','Trained','MM','Intermediate';
%     'CSP','Trained','MM','Deep'};

% toLoad = {'Plex','Trained','ALM','Superficial';...
%     'Plex','Trained','ALM','Intermediate';...
%     'Plex','Trained','ALM','Deep'};

% toLoad = {'Plex','Trained','MM','Superficial';...
%     'Plex','Trained','MM','Intermediate';...
%     'Plex','Trained','MM','Deep'};

toLoad = {'Fez','Trained','ALM','Deep';
    'Fez','Trained','MM','Deep'};

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
clear nVc aVc allVc idxRed;

S = twoP_settings;
sPerFrame = S.msPerFrame/1000;
trialType = 1; %1: stim, 2: response; 3: rewarded; 4: error

fData = length(E); 
idx = find(~cellfun(@isempty, E{1}(:,4)));
trialTypeNames = fieldnames(E{1}{idx(1),5}.sub.names);
fldBhv = fieldnames(E{1}{idx(1),5}.(trialTypeNames{trialType}));
posnVc = []; negNeurons = []; diffNeurons = [];
posallVc = []; negallVc = [];
posnVc = []; negnVc = [];
% meanVc = nan(length(fldBhv),)
for f = 1:length(E) % Number of cell type/cortical region/depth permutations for comparison
    allData = E{f};
    idxRow = find(~cellfun(@isempty, allData(:,4)));
    for j = 1:length(fldBhv)
        nVc = []; allVc = []; 
        idxRed{f} = [];
        countNan{f} = [];
        for i = 1:length(idxRow)
            tempnVc = mean(allData{idxRow(i),4}(:,:,allData{idxRow(i),5}.(trialTypeNames{trialType}).(fldBhv{j})),3,'omitnan');
            tempnVc(tempnVc > 3) = NaN;
            nVc = [nVc ; tempnVc];
            tempVc = mean(allData{idxRow(i),3}(:,:,allData{idxRow(i),5}.(trialTypeNames{trialType}).(fldBhv{j})),3,'omitnan');
            allVc = [allVc ; tempVc];
            idxRed{f} = [idxRed{f};allData{idxRow(i),6}];
            countNan{f} = [countNan{f};allData{idxRow(i),7}];
        end
        cellCount(f).notRed = sum(~idxRed{f}); cellCount(f).Red = sum(idxRed{f});
        allNan{f} = sum(countNan{f});
        filterNan = ones(1,length(allNan{f}));
        [valPks,idxPks] = findpeaks(allNan{f});
        mu = 0.5;
        for i = 1:length(idxPks)
            [valMin(i),idxMin(i)]=min(abs(allNan{f}(S.segFrames(i):idxPks(i))-mu*valPks(i)));
            filterNan(S.segFrames(i)+idxMin(i):idxPks(i))=nan;
        end
        t = 1:size(nVc,2);
        posnVc(j,:,f) = mean(nVc(find(idxRed{f}),:),1,'omitnan') .* filterNan;
        negnVc(j,:,f) = mean(nVc(find(~idxRed{f}),:),1,'omitnan') .* filterNan;
        posallVc(j,:,f) = mean(allVc(find(idxRed{f}),:),1,'omitnan') .* filterNan;
        negallVc(j,:,f) = mean(allVc(find(~idxRed{f}),:),1,'omitnan') .* filterNan;
        allVc = allVc.*filterNan;
        diffNeurons(j,:,f) = posnVc(j,:,f)-negnVc(j,:,f);
        aVc{f}(:,:,j) = allVc;
    end
end
meanVc = nan(size(posallVc,1), size(posallVc,2), size(posallVc,3),2);
meanVc(:,:,:,1) = posallVc;
meanVc(:,:,:,2) = negallVc;

%% Plot PETH of Vc -----------------------
close all;
t = 0:sPerFrame:sPerFrame*(size(meanVc,2)-1);
transparency = 0.1;
clrOrange = [0.9290 0.6940 0.1250]; clrMagenta = [1 0 1];
if size(meanVc,1) == 3
    lineColor = [[0 0 0]; clrMagenta; clrOrange];
    lineMarker = {'-';'-';'-'};
elseif size(meanVc,1) == 4
    lineColor = [clrMagenta; clrOrange; clrMagenta; clrOrange];
    lineMarker = {':';':';'-';'-'};
end

figMeanVc= figure(3);
set(figMeanVc,'Position',[100 100 400*size(meanVc,4) 300*size(meanVc,3)]);
tMeanVc = tiledlayout(size(meanVc,3),size(meanVc,4), ...
    'TileIndexing','ColumnMajor',...
    'TileSpacing','tight',...
    'Padding','tight');

% hold on;
% tt = 1;
bhvTrialTypeStart = size(aVc{1},3)-1;
for k = 1:2 % k refers to red versus non-red cells
%     YLim = squeeze([min(meanVc(:,:,:,k),[],[1 2],'omitnan') max(meanVc(:,:,:,k),[],[1 2],'omitnan')])';
    for i = 1:length(aVc) % i refers to the permutations for comparison
        lgdLines = [];
        nexttile
        for j = bhvTrialTypeStart:size(aVc{i},3)  % types of trials to be compared (e.g. right response)
            if k ==1
                [lineOut{j},~] = stdshade(aVc{i}(find(idxRed{i}),:,j),transparency,lineColor(j,:),t,[]);
            elseif k == 2
                [lineOut{j},~] = stdshade(aVc{i}(find(~idxRed{i}),:,j),transparency,lineColor(j,:),t,[]);                
            end
            %             line(t,meanVc(j,:,i,k),...
%                 'LineStyle',lineMarker{j},...
%                 'Color',lineColor(j,:),...
%                 'LineWidth',2);
            ax = gca;
            ax = fig_configAxis(ax);
            if i ~= length(aVc); ax.XTickLabel = []; ax.XAxis.Color = 'none'; end
            lgdLines = [lgdLines lineOut{j}];
        end
%         tt = tt+1;
        xlim([min(t) max(t)]);
        title({[M{i}.cellType ' ' M{i}.location ' ' M{i}.depth ' ' trialTypeNames{trialType} ' trials']});
        if k == 1 && i == 1
            title({[M{i}.cellType ' ' M{i}.location ' ' M{i}.depth ' ' trialTypeNames{trialType} ' trials'];'Labeled'});
        elseif k == 2 && i == 1
            title({[M{i}.cellType ' ' M{i}.location ' ' M{i}.depth ' ' trialTypeNames{trialType} ' trials'];'Unlabeled'});
        end
        ylabel({'Inferred spikes'});
        if i == 2 && k == 1
            lgd = legend(lgdLines, allData{idxRow(i),5}.sub.names.(trialTypeNames{trialType})(bhvTrialTypeStart:end),...
                'Location','Northwest',...
                'Box','off');
        end
        if i == 3; xlabel('Time (s)'); end
    end
end
exportgraphics(figMeanVc,fullfile(S.dir.imagingRootDir,[fName '_' trialTypeNames{trialType} '_IS.pdf']));

%% Plot epoch aligned Z-score_diff
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
set(figEpoch,'Position',[100 100 400 300*size(epochAligned,3)]);
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


%% Plot and save figure(s)

figAllPlots = figure(1);
set(figAllPlots,'Position',[100 100 1500 800]);
tFig = tiledlayout(length(E),length(fldBhv),...
    'TileSpacing','compact',...
    'Padding','compact');
leftYLim = squeeze([min([posnVc;negnVc],[],[1 2],'omitnan') max([posnVc;negnVc],[],[1 2],'omitnan')])';
rightYLim = squeeze([min(diffNeurons,[],[1 2],'omitnan') max(diffNeurons,[],[1 2],'omitnan')])';
t = 0:sPerFrame:sPerFrame*(size(posnVc,2)-1);

for f = 1:length(E)
    for j = 1:size(posnVc,1)
        nexttile
        
        hold on;
        
        yyaxis left
        line(t,posnVc(j,:,f),...
            'Color','r',...
            'LineWidth',1);
        line(t,negnVc(j,:,f),...
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