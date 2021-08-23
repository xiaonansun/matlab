function rateDisc_linePlot(ax,cData, cDates, cIdx)

plot(ax, cData(:, cDates >= cIdx(1) & cDates <= cIdx(2)),'linewidth',2,'color',[0.5 0.5 0.5]); hold on;
plot(ax, nanmean(cData(:, cDates >= cIdx(1) & cDates <= cIdx(2)),2),'linewidth',4,'color','k');