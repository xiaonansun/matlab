dataOverview = rateDiscRecordings;
animals = dataOverview(1:10,1);
groups = {'mSM', 'Fez', 'Plex' 'CSP'};
groupColor = {[1 0 0] [0 0 1] [0 0 0] [0 1 0]}; %colors for groups
cPath = 'R:\Behavior_Simon\'; %data path on the server

%% performance plot
% get nr of animals per group
clear groupCnt
for x = 1:length(groups)
    groupCnt{x} =find(contains(animals,groups{x}));
end

for iGroups = 1 : length(groups)
    figure('Renderer','painters');
    for iAnimals = 1 : length(groupCnt{iGroups})
        
        [Performance,bhv,~] = RateDisc_learningCurves(animals{groupCnt{iGroups}(iAnimals)},cPath,[],[],[],false);
        
        % learning curve plot
        subplot(1,3,1);
        cData = Performance.audioLearn.vals;
        cFit = Performance.audioLearn.cFit;
        fitRange = Performance.audioLearn.fitRange;
        cDate = Performance.audioLearn.dates;
        cColor = groupColor{iGroups} .* (((iAnimals / length(groupCnt{iGroups}))/2) + 0.5);

        plot(cDate, cData, '-', 'Color', cColor, 'MarkerFaceColor','w', 'linewidth', 1, 'MarkerSize', 4); hold on; ax = gca;
        plot(fitRange, cFit, 'Color', cColor, 'linewidth', 4); 
        ax.YLim = [0.4 1]; axis square; 
        xlabel('Sessions'); ylabel('Hit rate'); drawnow
        
        % discrimination over days
        subplot(1,3,2);
        cData = Performance.audioDisc.vals;
        cDate = Performance.audioDisc.dates;
        cDate = removeGaps(cDate(~isnan(cData)));
        cData = cData(~isnan(cData));

        a = unique(cDate);
        meanData = NaN(1, length(a));
        for x = 1:length(a)
            meanData(x) = mean(cData(cDate == a(x)));
        end
        plot(cDate, cData, '-', 'Color', cColor, 'MarkerFaceColor','w', 'linewidth', 1, 'MarkerSize', 4); hold on; ax = gca;
        plot(a, smooth(meanData,10,'lowess'), 'Color', cColor, 'linewidth', 4);
        ax.YLim = [0.4 1]; axis square;
        xlabel('Sessions'); ylabel('Hit rate'); drawnow
        
        % discrimination performance
        subplot(1,3,3);
        pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %only use active trials
        [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, pInd, 10);
        [~, ~, ~, ~, cFit] = rateDisc_audioDiscCurve(bhv, pInd, 10, true);

        %Plot fit and real data
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColor); hold on
        errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColor, 'MarkerFaceColor','w','linewidth',2);
        xlabel('right clicks'); ylabel('chose right'); axis square;
        ax = gca; ax.TickLength = [0 0];
        title('discrimination performance'); drawnow
        
    end
end
