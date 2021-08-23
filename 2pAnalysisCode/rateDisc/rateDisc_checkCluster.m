function rateDisc_checkCluster(T,X,mask,nrRows,savePath)
% code to find anatomical coordinates in the allen framework.
% reports coordinates in mm, relative to bregma location.

h = figure('Renderer','painters');
h.UserData.T = T;
h.UserData.lastT = T;
h.UserData.X = X;
h.UserData.mask = mask;
h.UserData.nrRows = nrRows;
h.UserData.savePath = savePath;
h.UserData.cClust = [];

h.UserData.corrFig = figure;
h.UserData.corrMat = corrcoef(h.UserData.X');
updateCorrFig(h);

%% initial plots
drawClusters(h);

% merge button
ax = h.Children(end);
uicontrol('Style', 'pushbutton', 'String', 'Merge', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) + 0.1, 0.1 0.04], ...
    'Callback', @clusterMerge, 'Enable', 'on');

% partial merge button
ax = h.Children(end);
uicontrol('Style', 'pushbutton', 'String', 'partialMerge', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) + 0.05, 0.1 0.04], ...
    'Callback', @clusterPartMerge, 'Enable', 'on');

% split button
ax = h.Children(end);
uicontrol('Style', 'pushbutton', 'String', 'split', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) + 0, 0.1 0.04], ...
    'Callback', @clusterSplit, 'Enable', 'on');

% undo button
ax = h.Children(end);
uicontrol('Style', 'pushbutton', 'String', 'undo merge/split', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) - 0.05, 0.1 0.04], ...
    'Callback', @undoMerge, 'Enable', 'on');

% check button
ax = h.Children(end);
uicontrol('Style', 'pushbutton', 'String', 'show members', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) - 0.1, 0.1 0.04], ...
    'Callback', @showMembers, 'Enable', 'on');

% reset button
uicontrol('Style', 'pushbutton', 'String', 'Reset', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) - 0.15, 0.1 0.04], ...
    'Callback', @clusterReset, 'Enable', 'on');

% done button
uicontrol('Style', 'pushbutton', 'String', 'done', ...
    'Units', 'Normalized', 'Position', [ax.Position(1) - 0.11, ax.Position(2) - 0.2, 0.1 0.04], ...
    'Callback', @clusterOut, 'Enable', 'on');


% uiwait(h);
end

%% outline callback
function showOutline(hObject,~,~)
ax = hObject.Parent;
checker = textscan(ax.Title.String, '%f%s%f');
checker = checker{1};
ax.Parent.UserData.cClust = checker;

if ax.LineWidth ~= 5
    ax.LineWidth = 5;
    grid(ax,'on');  
    if strcmpi(ax.Parent.SelectionType, 'normal')
        ax.XColor = 'r'; ax.YColor = 'r';
    else
        ax.XColor = 'b'; ax.YColor = 'b';
    end
    
    % mark in correlation matrix
    [a,b] = sort(ax.Parent.UserData.T);
    cOn = find(a==checker,1,'first');
    cSize = find(a==checker,1,'last') - cOn + 1;
    
    figure(ax.Parent.UserData.corrFig);
    axes(ax.Parent.UserData.corrMatAx);
    ax.Parent.UserData.markers(checker) = rectangle('Position',[cOn, cOn, cSize, cSize], 'EdgeColor', 'k', 'LineWidth',4);
    
    % show correlation with other clusters
    cla(ax.Parent.UserData.clustcorrAx);
    corrMat = ax.Parent.UserData.corrMat(b,b);
    cData = corrMat(cOn : cOn + cSize-1, :)';
    cData = cData(:);
    axes(ax.Parent.UserData.clustcorrAx); hold on;
    plot([find(unique(a)==checker) find(unique(a)==checker)],[0 1.1], 'k');
    plot([0 length(unique(a))+1],[0.9 0.9], 'k');   
    violinplot(cData,repmat(a,cSize,1),'ViolinAlpha',0.25,'ShowData', false, 'ShowMean', true);
    xlabel('Clusters'); ylabel('Correlation'); title('Correlation with other clusters'); 
    ax.Parent.UserData.clustcorrAx.TickLength = [0 0];
    grid on;

else
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 0.5;
    delete(ax.Parent.UserData.markers(checker));
    cla(ax.Parent.UserData.clustcorrAx);
end
end

%% reset button callback
function clusterReset(hObject,~,~)
ax = hObject.Parent;
delete(ax.UserData.markers(:));
cla(ax.UserData.clustcorrAx);
for x = 1 : length(ax.Children)
    try
        if ax.Children(x).LineWidth == 5
            ax.Children(x).LineWidth = 0.5;
            ax.Children(x).XColor = 'k'; ax.Children(x).YColor = 'k';
        end
    end
end
end

%% merge button callback
function clusterMerge(hObject,~,~)
ax = hObject.Parent;
cSelect = [];
ax.UserData.lastT = ax.UserData.T;

for x = 1 : length(ax.Children)
    checker = 0;
    try
        checker = textscan(ax.Children(x).Title.String, '%f%s%f');
        checker = checker{1};
    end
    if checker > 0
        if ax.Children(x).LineWidth == 5
            cSelect = [cSelect checker];
            ax.Children(x).LineWidth = 0.5;
            ax.Children(x).XColor = 'k'; ax.Children(x).YColor = 'k';
        end
    end
end

if ~isempty(cSelect)
    % merge clusters
    ax.UserData.T(ismember(ax.UserData.T,cSelect)) = min(cSelect);
    updateCorrFig(ax);
    drawClusters(ax); %redraw clusters
end
end

%% partial merge button callback
function clusterPartMerge(hObject,~,~)
ax = hObject.Parent;
cSelect = [];
ax.UserData.lastT = ax.UserData.T;

for x = 1 : length(ax.Children)
    checker = 0;
    try
        checker = textscan(ax.Children(x).Title.String, '%f%s%f');
        checker = checker{1};
    end
    if checker > 0
        if ax.Children(x).LineWidth == 5
            if ax.Children(x).XColor(3) == 1
                cSelect = [checker cSelect]; %make this the main cluster (by default its the lowest number)
            else
                cSelect = [cSelect checker];
            end
            ax.Children(x).LineWidth = 0.5;
            ax.Children(x).XColor = 'k'; ax.Children(x).YColor = 'k';
        end
    end
end

if ~isempty(cSelect)
    %if only one cluster is selected, search for similar dims in other clusters
    if length(cSelect) == 1 
        cIdx = true(1,length(ax.UserData.T));
    else %merge only selected clusters
        cIdx = ismember(ax.UserData.T,cSelect);
    end
    
    % merge clusters
    meanCluster = nanmean(ax.UserData.X(ax.UserData.T == cSelect(1),:),1);
    xx = corrcoef([meanCluster; ax.UserData.X(cIdx,:)]');
    if length(cSelect) == 1 
        cIdx(cIdx) = xx(1,2:end) > 0.9; %change index to remove components below threshold
    else
        cIdx(cIdx) = xx(1,2:end) > 0.9; %for partial merge between selected clusters, use lower threshold
    end
    ax.UserData.T(cIdx) = cSelect(1);
    updateCorrFig(ax);
    drawClusters(ax); %redraw clusters
end
end


%% split button callback
function clusterSplit(hObject,~,~)
ax = hObject.Parent;
cSelect = [];
ax.UserData.lastT = ax.UserData.T;

for x = 1 : length(ax.Children)
    checker = 0;
    try
        checker = textscan(ax.Children(x).Title.String, '%f%s%f');
        checker = checker{1};
    end
    if checker > 0
        if ax.Children(x).LineWidth == 5
            cSelect = [cSelect checker];
            ax.Children(x).LineWidth = 0.5;
            ax.Children(x).XColor = 'k'; ax.Children(x).YColor = 'k';
        end
    end
end

if length(cSelect) > 1
    disp('More than one cluster select. Can only split one!');
elseif isempty(cSelect)
    disp('No cluster selected!');
else
    % split cluster
    cIdx = ax.UserData.T == cSelect(1);
    newClust = kmeans(ax.UserData.X(cIdx,:),2)+1000;
    temp = unique(ax.UserData.T)';
    if any(ismember(temp,cSelect(1)+1)) %make space for new cluster if needed
        ax.UserData.T(ax.UserData.T > cSelect(1)) = ax.UserData.T(ax.UserData.T > cSelect(1)) + 1;
    end
    newClust(newClust == 1001) = cSelect(1); %old cluster
    newClust(newClust == 1002) = cSelect(1)+1; %new cluster
    ax.UserData.T(cIdx) = newClust;
    updateCorrFig(ax);
    drawClusters(ax); %redraw clusters
end
end

%% undo button
function undoMerge(hObject,~,~)
    % merge clusters
    ax = hObject.Parent;
    ax.UserData.T = ax.UserData.lastT;
    updateCorrFig(ax);
    drawClusters(ax); %redraw clusters
end

%% show members
function showMembers(hObject,~,~)
    if ~isempty(hObject.Parent.UserData.cClust)
        cIdx = hObject.Parent.UserData.T == hObject.Parent.UserData.cClust; %members of current cluster
        figure; imagesc(corrcoef(hObject.Parent.UserData.X(cIdx,:)')); axis square; colormap viridis; colorbar
        compareMovie(arrayShrink(hObject.Parent.UserData.X(cIdx,:)',hObject.Parent.UserData.mask,'split'));
    end
end

%% done button callback
function clusterOut(hObject,~,~)
ax = hObject.Parent;
clusters = cell(1,length(unique(ax.UserData.T)));
for x = 1 : length(ax.Children)
    if contains(class(ax.Children(x)),'Axes')
        try
            checker = textscan(ax.Children(x).Title.String, '%f%s%f');
            clusters{checker{1}} = ax.Children(x).Children.CData;
            if ax.Children(x).LineWidth == 5
                ax.Children(x).LineWidth = 0.5;
                ax.Children(x).XColor = 'k'; ax.Children(x).YColor = 'k';
            end
        end
    end
end
% uiresume(ax);
newT = ax.UserData.T;
assignin('base','newT',newT); %send current image of both axes in cell array. Cell 1 is leftaxis , cell 2 is right axis.
corrMat = ax.UserData.corrMat; %keep sorted correlation matrix
assignin('base','corrMat',corrMat); %send correlation matrix to workspace

cFile = ax.UserData.savePath;
if ~exist(fileparts(cFile),'dir')
    mkdir(fileparts(cFile));
end
if exists(cFile) %if file exists already keep a copy of the original
    mkdir([fileparts(cFile) filesep 'old']);
    movefile(cFile, strrep(cFile, fileparts(cFile), [fileparts(cFile) filesep 'old']));
end
save(cFile, 'newT', 'clusters','corrMat');

% close(ax.UserData.corrFig);
% close(ax);
end


%% initial plots
function drawClusters(h)

figure(h);
cClust = unique(h.UserData.T)';
rejIdx = [];
for x = 1 : length(h.Children)
    if contains(class(h.Children(x)),'Axes')
        rejIdx = [rejIdx x];
    end
end
delete(h.Children(rejIdx)); %remove old plots

Cnt = 0;
for x = cClust
    Cnt = Cnt + 1;
    cIdx = h.UserData.T == x;
    subplot(h.UserData.nrRows, ceil(length(cClust)/h.UserData.nrRows),Cnt);
    cImg = imagesc(arrayShrink(nanmean(h.UserData.X(cIdx,:),1)',h.UserData.mask,'split'),'ButtonDownFcn',{@showOutline}); axis image; colormap colorcube
    ax = cImg.Parent; grid(ax,'on');ax.GridColor = 'w';
    ax.XTickLabel = []; ax.YTickLabel = [];
    ax.BoxStyle = 'full'; ax.Box = 'on';
    title([num2str(x) ' - corr: ' num2str(h.UserData.mCorr(Cnt)) '; Cnt: ' num2str(sum(cIdx))]); axis image; colormap('jet');
end
end

%% update correlation matrix
function updateCorrFig(h)

figure(h.UserData.corrFig); hold off;
[a,b] = sort(h.UserData.T);
xx = h.UserData.corrMat(b,b);
h.UserData.corrMatAx = subplot(6,2,1:10);
imagesc(smoothImg(xx)); axis square
colormap viridis; 
caxis([0.5 0.9]);
h.UserData.mCorr = NaN(1, length(unique(h.UserData.T)));

Cnt = 0;
for x = unique(h.UserData.T)'
    try
        cIdx = a==x;
        vline(find(cIdx,1,'last'),'w'); text(find(cIdx,1,'last')-10, 30, num2str(x), 'Color', 'w')
        hline(find(cIdx,1,'last'),'w'); text(30, find(cIdx,1,'last')-10, num2str(x), 'Color', 'w')
        Cnt = Cnt + 1;
        h.UserData.mCorr(Cnt) = round(nanmean(reshape(xx(cIdx,cIdx),1,[])),2);
    end
end

% to indicate correlation with other clusters
h.UserData.clustcorrAx = subplot(6,2,11:12);

end