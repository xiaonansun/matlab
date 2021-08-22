function twoP_plotAUCFraction(A,AUC)

figureTitle = [A.currentAnimal ' ' AUC.location];
figureSaveFilename = fullfile(A.Setting.baseDir,[A.currentAnimal '_' AUC.location '_' A.allAUC.analysisType '_allDistributions.pdf']);
P(1).traceFraction = AUC.R.frPos; P(1).mean = AUC.R.meanPos; P(1).legend = 'Ipsi Red'; P(1).color = 'r'; P(1).linestyle = '-';
P(2).traceFraction = AUC.NR.frPos; P(2).mean = AUC.NR.meanPos; P(2).legend = 'Ipsi non-Red'; P(2).color = 'k'; P(2).linestyle = '-';
P(3).traceFraction = AUC.R.frNeg; P(3).mean = AUC.R.meanNeg; P(3).legend = 'Contra Red'; P(3).color = 'r'; P(3).linestyle = ':';
P(4).traceFraction = AUC.NR.frNeg; P(4).mean = AUC.NR.meanNeg; P(4).legend = 'Contra non-Red'; P(4).color = 'k'; P(4).linestyle = ':';
segFrames = [30    53    91   106   136]; % this is temporary, need to find a way to put this into the struct so it can be loaded
frames = 1:1:length(AUC.R.frPos);

figTimeCourseAUC = figure(3);
set(figTimeCourseAUC, 'Position', [100 100 500 500]);

t=tiledlayout(2,1);
nexttile; hold on;
for i = 1:length(P)
    line(frames, P(i).traceFraction, ...
    'LineWidth',2,...
    'LineStyle',P(i).linestyle,...
    'Color',P(i).color);
end
xline(segFrames+1);
ax = gca;
ax = fig_configAxis(ax);
title(figureTitle);  ylabel('Fraction choice selective')
legend({P(:).legend},'Location','northwest');

nexttile; hold on;
for i = 1:length(P)
    line(frames, P(i).mean, ...
    'LineWidth',2,...
    'LineStyle',P(i).linestyle,...
    'Color',P(i).color);
end
ylim([0.2 0.8]);
xline(segFrames+1);
ax = gca;
ax = fig_configAxis(ax);
xlabel('Frame'); ylabel('Mean AUC of choice-selective neurons')

exportgraphics(gcf,figureSaveFilename);