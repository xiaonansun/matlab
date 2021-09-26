if ~exist('Vsub','var')
    load('\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy\analysis\all_psth.mat','Vsub','trialTypes');
end
clear fPETH fDiff;

S = twoP_settings;
y=struct; d = struct;
P = cell(1,length(S.cellTypes));
SP = cell(1,length(S.cellTypes));
catCR = cell(1,length(S.cellTypes));
catCNR = cell(1,length(S.cellTypes));
catSR = cell(1,length(S.cellTypes));
catSNR = cell(1,length(S.cellTypes));

% trialTypes = {'StimAllRight','StimAllLeft','StimEasyRight','StimEasyLeft','ResponseAll','ResponseLeft','ResponseRight','RewardedAll','RewardedLeft','RewardedRight','ErrorAll','ErrorLeft','ErrorRight'};
lineWidth = 1;

noCharIdx = cellfun(@isempty,Vsub(:,end));
Vsub{noCharIdx,end} = '';

for i = 1:length(P)
    clear ctIdx;
    ctIdx = contains(Vsub(:,end),S.cellTypes{i}); %% ctIdx: cell type index
    ctIdx = ~cellfun(@isempty,Vsub(:,1)) & ctIdx;
    CR = Vsub(ctIdx,1); CNR= Vsub(ctIdx,3); 
    SR = Vsub(ctIdx,5); SNR = Vsub(ctIdx,7);
    catCR{i} = vertcat(CR{:}); catCNR{i} = vertcat(CNR{:});
    catSR{i} = vertcat(SR{:}); catSNR{i} = vertcat(SNR{:});
    P{i} = vertcat(catCR{i},catCNR{i});
    SP{i} = vertcat(catSR{i},catSNR{i});
    y(i).R = squeeze(mean(catCR{i},1,'omitnan')); y(i).NR = squeeze(mean(catCNR{i},1,'omitnan'));
    y(i).SR = squeeze(mean(catSR{i},1,'omitnan')); y(i).SNR = squeeze(mean(catSNR{i},1,'omitnan'));
    y(i).P = squeeze(mean(P{i},1,'omitnan')); y(i).SP = squeeze(mean(SP{i},1,'omitnan'));
    d(i).R = y(i).R-y(i).SR; d(i).NR = y(i).NR-y(i).SNR;
    d(i).P = y(i).P - y(i).SP;
    temp = struct2cell(y(i)); temp = vertcat(temp{:});
    yLimit(i,:)= [min(temp(:)) max(temp(:))]; clear temp
    temp = struct2cell(d(i)); temp = vertcat(temp{:});
    dLimit(i,:)= [min(temp(:)) max(temp(:))]; clear temp
end
x = 1:1:size(P{1},2);


fPETH = figure(1);
set(fPETH, 'Position',[0 0 2000 1000]);
hT = tiledlayout(length(P),size(P{1},3),...
    'TileSpacing','tight',...
    'Padding','tight');
for iF = 1:length(P)
    for iT = 1:size(P{iF},3)
        nexttile; hold on;
        line(x,y(iF).R(:,iT),'Color','r',...
            'LineWidth',lineWidth); 
        line(x,y(iF).NR(:,iT),'Color','g',...
            'LineWidth',lineWidth);
        line(x,y(iF).SR(:,iT),...
            'Color','r',...
            'LineStyle','--',...
            'LineWidth',lineWidth);
        line(x,y(iF).SR(:,iT),...
            'Color','g',...
            'LineStyle','--',...
            'LineWidth',lineWidth);
        line(x,y(iF).P(:,iT),...
            'Color','k',...
            'LineStyle','-',...
            'LineWidth',lineWidth);
        line(x,y(iF).SP(:,iT),...
            'Color','k',...
            'LineStyle','--',...
            'LineWidth',lineWidth);
        ylim([yLimit(iF,1) yLimit(iF,2)]);
        ax = gca;
        ax = fig_configAxis(ax);
        if iF == 1; title(trialTypes{2,iT}); end
        if iF ~= 3; ax.XTickLabel = []; ax.XAxis.Color = 'none'; end
        if iT == 1; ylabel(S.cellTypes{iF}); end
        if iT ~= 1; ax.YTickLabel = []; ax.YAxis.Color = 'none'; end
    end
end
xlabel(hT,'Frames'); ylabel(hT,'Inferred spikes');
title(hT,'PETH');
exportgraphics(fPETH,fullfile(S.dir.imagingRootDir,'analysis','PETH_by_trial_type.pdf'));


fDiff = figure(2);
set(fDiff, 'Position',[0 0 2000 1000]);
hTdiff = tiledlayout(length(P),size(P{1},3),...
    'TileSpacing','tight',...
    'Padding','tight');
for iF = 1:length(P)
    for iT = 1:size(P{iF},3)
        nexttile; hold on;
        line(x,d(iF).R(:,iT),'Color','r',...
            'LineWidth',lineWidth);
        line(x,d(iF).NR(:,iT),'Color','g',...
            'LineWidth',lineWidth);
        line(x,d(iF).P(:,iT),'Color','k',...
            'LineWidth',lineWidth);
        yline(0);
        ylim([dLimit(iF,1) dLimit(iF,2)]);
        ax = gca;
        ax = fig_configAxis(ax);
        if iF == 1; title(trialTypes{2,iT}); end
        if iF ~= 3; ax.XTickLabel = []; ax.XAxis.Color = 'none'; end
        if iT == 1; ylabel(S.cellTypes{iF}); end
        if iT ~= 1; ax.YTickLabel = []; ax.YAxis.Color = 'none'; end
    end
end
xlabel(hTdiff,'Frames'); ylabel(hTdiff,'Inferred spikes');
title(hTdiff,'PETH - shuffled PETH');
exportgraphics(fDiff,fullfile(S.dir.imagingRootDir,'analysis','PETH_by_trial_type_minus_shuffle.pdf'));