cPath = '\\grid-hs\churchland_nlsas_data\BehaviorVideo\';
animal = 'CSP29';
fPath = [cPath animal filesep 'SpatialDisc' filesep 'Session Data' filesep];
earlyDate = 'Feb29_2020';

%% load behavior files
bhvFiles = dir([fPath animal '*.mat']);
bhv = [];

for iFiles = 1 : length(bhvFiles)
    SessionData = [];
    load([fPath bhvFiles(iFiles).name], 'SessionData');
    
    if ~isempty(SessionData)
        bhv = appendBehavior(bhv,SessionData); %append into larger array
    end
end

%% compute discrimination performance.
distBins = 6;
cIdx(1,:) = bhv.TrialStartTime < datenum(earlyDate);
cIdx(2,:) = bhv.TrialStartTime > datenum(earlyDate);
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %only use active trials

figure
cColors = {'k', 'r'};
cTitles ={'naive', 'expert'};
clear h
for x = 1 : 2
    
    subplot(1,2,x);
    [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, cIdx(x,:) & pInd, distBins);
    
    %Plot fit and real data
    plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
    h(x) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
    
    % add labels
    xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
    ylim([0 1]); ylabel('Proportion chose right'); hold off; axis square
    xlabel('Distractor ratio'); title([animal ' - ' cTitles{x}]);
    disp(nTrials);
    
end