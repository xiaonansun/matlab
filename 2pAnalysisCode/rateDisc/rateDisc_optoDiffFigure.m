function rateDisc_optoDiffFigure(bhv,groupnames)
%code to produce the panels about performance change over the course of
%training. Compares impairment in detection performance early and late in
%training for different cell types. expects bhv to be a cell array where
%every cell is a cell-type.

h = figure; 
fiberColors = {[0 0 1] [1 0 0]};
Cnt = 0;
for x = 1 : length(bhv)
    for cLoc = 1 : 2 %location (1 = frontal, 2 = parietal)
        clear cData
        out = rateDisc_stimANDdelayOptoDiff(bhv{x}); %out.detect: early/late x ctrl/opto x location x animals
        
        ax = gca; hold on;
        cData = squeeze(out.detect(:,1,cLoc,:) - out.detect(:,2,cLoc,:)).*100; %difference between control and stimulation trials
        cData = cData(:,~isnan(mean(cData))); %only use mice that have early and late recordings
        errorbar(Cnt : Cnt+1, nanmean(cData,2), nansem(cData,2), 'k.', 'linewidth', 2)
        bar(Cnt : Cnt+1,nanmean(cData,2), 'FaceColor', fiberColors{cLoc}, 'FaceAlpha', 0.5)
        plot(Cnt : Cnt+1, cData, '-o', 'Color', [fiberColors{cLoc} 0.25], 'LineWidth', 2, 'MarkerFaceColor', 'w')        
        
%         boxplot(squeeze(out.detect(1,1,cLoc,:) - out.detect(1,2,cLoc,:))'.*100, 'positions', Cnt, 'Colors', fiberColors{cLoc})
%         boxplot(squeeze(out.detect(2,1,cLoc,:) - out.detect(2,2,cLoc,:))'.*100, 'positions', Cnt+1, 'Colors', fiberColors{cLoc})
        
%         Violin(cData(1,:)', Cnt, 'ShowData', false, ...
%             'ViolinAlpha',0.1, 'ViolinColor', fiberColors{cLoc}, 'Width', 0.2, 'Bandwidth', 5, 'MedianColor', fiberColors{cLoc}); %early
%         
%         Violin(cData(1,:)', Cnt+1, 'ShowData', false, ...
%             'ViolinAlpha',0.1, 'ViolinColor', fiberColors{cLoc}, 'Width', 0.2, 'Bandwidth', 5, 'MedianColor', fiberColors{cLoc}); %late
        Cnt = Cnt + 2;
    end
    Cnt = Cnt + 2;
end
ax.XTick = 1.5:6:Cnt;
ax.XTickLabel = groupnames;
ax.TickLength = [0 0];
nhline(0, '--', 'Color', [0.5 0.5 0.5]); grid on;
ylabel('Detection impairment (%)');
xlim([-1 Cnt-1]); axis square; title('Performance change - Early vs late sessions');