function rateDisc_optoBasicFigure(bhv,aLabel)
%code to to show the general impact of optgenetic inactivation during the
%stimulus and delay period for frontal and parietal cortex. bhv is a cell
%array with behavioral data, each entry is a different cell-type.

%% compute basic inactivation effect during stimulus+delay period, early in training
h = figure('name', aLabel);
out = rateDisc_stimANDdelayOptoDiff(bhv); %out.detect: early/late x ctrl/opto x location x animals

subplot(1,2,1);
cData = squeeze(out.detect(1,1,:,:));
cData = [nanmean(cData(:)); nanmean(squeeze(out.detect(1,2,:,:)),2)];

b(1) = bar(0,cData(1),1,'facecolor',[0.5 0.5 0.5]); hold on
plot(0,nanmean(squeeze(out.detect(1,1,:,:)),1), 'ko','MarkerFaceColor','w');
b(2) = bar(1.5,cData(2),1);
plot(1.5,squeeze(out.detect(1,2,1,:)), 'ko','MarkerFaceColor','w');
b(3) = bar(3,cData(3),1);
plot(3,squeeze(out.detect(1,2,2,:)), 'ko','MarkerFaceColor','w');
xLabel = {'Control' 'Frontal' 'Parietal'};
grid on; b(1).Parent.XTick = [0 1.5 3]; b(1).Parent.XTickLabel = xLabel;
xlim([-1 4]); ylabel('Detection performance');
ylim([0.5 1]);  axis square; title(['Detection, Stim+Delay, ' aLabel]);
niceFigure(gca)


%% compute detection performance with bilateral during different task episodes from handle to response period
xLabels =  {'Handle' 'EarlyStim' 'LateStim' 'Delay' 'Response'};
fiberLocations = { 'Frontal' 'Parietal'};
fiberColors = {[0 0 1] [1 0 0]};
nrMice = length(unique(bhv.AnimalID));
optoPower = 10;

allData = rateDisc_taskEpisodesOpto(bhv, bhv.optoPower == optoPower); % get task episode data, low power
nGroups = size(allData.optoDetect,2);
allP = NaN(nGroups, size(allData.optoDetect,3), nrMice); %times x locations x animals

subplot(1,2,2); hold on
for x = 2 : -1 : 1
    for iAnimals = unique(bhv.AnimalID)
        allData = rateDisc_taskEpisodesOpto(bhv, bhv.optoPower == optoPower & ismember(bhv.AnimalID,iAnimals)); % get task episode data, low power
        dOut.allTimes{x,iAnimals} = allData; %keep for output
        allP(:,x,iAnimals) = allData.detect(3) - allData.optoDetect(3,:,x);
        plot(squeeze(allP(:,x,iAnimals)), 'Color', [fiberColors{x} 0.2]);
    end
    allData = rateDisc_taskEpisodesOpto(bhv, bhv.optoPower == optoPower); % get task episode data, low power
    dOut.allTimes{x,iAnimals + 1} = allData; %keep for output
    pChange = allData.detect(3) - [allData.optoDetect(3,:,x); allData.optoDetectUp(3,:,x); allData.optoDetectLow(3,:,x)];
    cLine(x) = errorbar(1:nGroups, pChange(1,:), pChange(1,:) - pChange(2,:), pChange(3,:) - pChange(1,:), '-o', 'linewidth' ,4, 'color', fiberColors{x},'MarkerFaceColor','w', 'MarkerSize', 10);
end

xlim([0.5 cLine(1).XData(end)+0.5]);
nhline(nanmean(0), '--', 'lineWidth',4, 'Color', [0.5 0.5 0.5]);
axis square;  
ylim([-0.1 0.3]);

cLine(1).Parent.XTick = 1:nGroups;
grid on;  cLine(1).Parent.XTickLabel = arrayfun(@(y) sprintf('%s - %i', xLabels{y}, allData.indCnt(y,x)), 1:5, 'UniformOutput',false);
title([aLabel ' (' num2str(optoPower) 'mW)']);
ylabel('Detection impairment (%)');
set(h,'PaperOrientation','landscape','PaperPositionMode','auto');
niceFigure(gca)
legend(cLine, fiberLocations);


%% add induction date for control mice
if any(ismember(bhv.Animals, 'Fez7'))
    minDate = 737549; % last date for Fez7 before being induced with TMX (first injection was on 5/4/19
else
    minDate = 1;
end

%% control figure for Fez7
if any(ismember(bhv.Animals, 'Fez7'))
    h = figure;
    %before induction
    allData = rateDisc_stimANDdelayOpto(bhv,bhv.date < minDate & ismember(bhv.AnimalID,find(ismember(bhv.Animals, 'Fez7')))); % get task episode data, low power
    subplot(1,2,1);
    errorbar(0, allData.detect(3), allData.detect(3) - allData.detectLow(3), allData.detectUp(3) - allData.detect(3),'-k', 'linewidth',2); hold on;
    h(1) = bar(0,allData.detect(3),1,'facecolor',[0.5 0.5 0.5]); hold on
    errorbar(1.5, allData.optoDetect(3,1), allData.optoDetect(3,1) - allData.optoDetectLow(3,1), allData.optoDetectUp(3,1) - allData.optoDetect(3,1),'-k', 'linewidth',2);
    h(2) = bar(1.5,allData.optoDetect(3,1),1);
    errorbar(3, allData.optoDetect(3,2), allData.optoDetect(3,2) - allData.optoDetectLow(3,2), allData.optoDetectUp(3,2) - allData.optoDetect(3,2),'-k', 'linewidth',2);
    h(3) = bar(3,allData.optoDetect(3,2),1);
    errorbar(4.5, allData.optoDetect(3,3), allData.optoDetect(3,3) - allData.optoDetectLow(3,3), allData.optoDetectUp(3,3) - allData.optoDetect(3,3),'-k', 'linewidth',2);
    h(4) = bar(4.5,allData.optoDetect(3,3),1);
    
    stimCnt = textscan(num2str(allData.optoCnt),'%s%s%s%s');
    xLabel = {['Control\newline' stimCnt{4}{1}] ['Frontal\newline' stimCnt{1}{1}] ['Parietal\newline' stimCnt{2}{1}]};
    grid on; h(1).Parent.XTick = [0 1.5 3]; h(1).Parent.XTickLabel = xLabel;
    xlim([-1 4]); ylabel('Detection performance');
    ylim([0.5 1]);  axis square; title('Detection, Stim+Delay inactivation - BEFORE TAMOXIFEN');    

    % after induction
    allData = rateDisc_stimANDdelayOpto(bhv,bhv.date > minDate & ismember(bhv.AnimalID,find(ismember(bhv.Animals, 'Fez7')))); % get task episode data, low power
    subplot(1,2,2);
    errorbar(0, allData.detect(3), allData.detect(3) - allData.detectLow(3), allData.detectUp(3) - allData.detect(3),'-k', 'linewidth',2); hold on;
    h(1) = bar(0,allData.detect(3),1,'facecolor',[0.5 0.5 0.5]); hold on
    errorbar(1.5, allData.optoDetect(3,1), allData.optoDetect(3,1) - allData.optoDetectLow(3,1), allData.optoDetectUp(3,1) - allData.optoDetect(3,1),'-k', 'linewidth',2);
    h(2) = bar(1.5,allData.optoDetect(3,1),1);
    errorbar(3, allData.optoDetect(3,2), allData.optoDetect(3,2) - allData.optoDetectLow(3,2), allData.optoDetectUp(3,2) - allData.optoDetect(3,2),'-k', 'linewidth',2);
    h(3) = bar(3,allData.optoDetect(3,2),1);
    errorbar(4.5, allData.optoDetect(3,3), allData.optoDetect(3,3) - allData.optoDetectLow(3,3), allData.optoDetectUp(3,3) - allData.optoDetect(3,3),'-k', 'linewidth',2);
    h(4) = bar(4.5,allData.optoDetect(3,3),1);
    
    stimCnt = textscan(num2str(allData.optoCnt),'%s%s%s%s');
    xLabel = {['Control\newline' stimCnt{4}{1}] ['Frontal\newline' stimCnt{1}{1}] ['Parietal\newline' stimCnt{2}{1}]};
    grid on; h(1).Parent.XTick = [0 1.5 3]; h(1).Parent.XTickLabel = xLabel;
    xlim([-1 4]); ylabel('Detection performance');
    ylim([0.5 1]);  axis square; title('Detection, Stim+Delay inactivation - AFTER TAMOXIFEN');    
end

