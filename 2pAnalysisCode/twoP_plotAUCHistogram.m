function twoP_plotAUCHistogram(A,AUC)

figureTitle = [A.currentAnimal ' ' AUC.location];
figureSaveFilename = fullfile(A.Setting.baseDir,[A.currentAnimal '_' AUC.location '_' A.allAUC.analysisType '_allDistributions.pdf']);
nBins = 20;
figAUC = figure(2);
set(figAUC, 'Position', [100 100 800 200]);
% hS = histogram(S.NR,20,'Normalization','probability');subplot(2,1,1);
v.NR.val = nan(size(AUC.NR.val,2),nBins); edges.NR.val= nan(size(AUC.NR.val,2),nBins+1);
for j = 1:size(AUC.NR.val,2)
[v.NR.val(j,:), edges.NR.val(j,:)] = histcounts(AUC.NR.val(:,j),nBins,'Normalization','probability');
c.NR.val(j,:) = conv(edges.NR.val(j,:),[0.5 0.5],'valid'); % computes the centers of each bin.
c.NR.view = [-90 90];
[v.R.val(j,:), edges.R.val(j,:)] = histcounts(AUC.R.val(:,j),nBins,'Normalization','probability'); 
c.R.val(j,:) = conv(edges.R.val(j,:),[0.5 0.5],'valid'); % computes the centers of each bin
c.R.view = [90 -90];
end

t=tiledlayout(1,2*(size(AUC.NR.val,2)));
m = 1; n = 1; maxYLim = 0.4;
for i = 1: 2*size(AUC.NR.val,2)
    if ~rem(i,2)==0
        nexttile; hold on;
        line(c.NR.val(m,:),v.NR.val(m,:),...
            'LineWidth',2,...
            'Color','k'); view(c.NR.view);
        plot([0.5 0.5],[0 0.5],...
            'color',[0.5 0.5 0.5]);
        title(AUC.epochNames{m}(2:end));
        ax = gca; 
        if m == 1; xlabel('AUC');  end
        if m > 1; ax.XAxis.Visible = 'off';  end
        ax = fig_configAxis(ax);
        xlim([0 1]); ylim([0 maxYLim]);
        m = m+1;
    elseif rem(i,2)==0
        nexttile; hold on;
        line(c.R.val(n,:), v.R.val(n,:),...
            'LineWidth',2,...
            'Color','g'); view(c.R.view);
        plot([0.5 0.5],[0 0.5],...
            'color',[0.5 0.5 0.5]);
        ax=gca; ax.XAxis.Visible = 'off';
        if n >= 1; end
        ax = fig_configAxis(ax);
        xlim([0 1]), ylim([0 maxYLim])
        n = n+1;
    end
end

t.Padding = 'tight'; t.TileSpacing = 'none';
title(t,figureTitle); xlabel(t, 'p')
exportgraphics(gcf,figureSaveFilename);
