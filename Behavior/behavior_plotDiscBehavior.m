%%
close all;
Animal = 'Plex62';
cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
[Performance,bhv] = DelayedLoc_learningCurves(Animal,cPath);
disp([Animal '. Detection sessions: ' num2str(sum(sum(Performance.Detection>0))) '; Discrimination Sessions: ' num2str(sum(sum(Performance.Discrimination>0)))]);
% [Performance,bhv,allDates] = RateDisc_learningCurves(Animal,cPath);

%% Plot discrimination curve
hDisc = figure(3);
cInd = 1:length(SessionData.Rewarded);
[distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(SessionData, cInd);
% rateDisc_audioDiscCurve(bhv, cInd, distBins, discOnly, fixBias, returnCIs)
plot(distRatio,pChoseHigh,'.-k'); hold on;
plot(distRatio,dataUpper,'.--k'); 
plot(distRatio,dataLower,'.--k'); 

xlim([0 1]); ylim([0 1]);
title([num2str(sum(SessionData.Rewarded)) ' trials rewarded. ' num2str(length(SessionData.Rewarded)) ' trials completed. ']);
ylabel('Proportion of right choices');

%%
cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
animal = 'Plex62';

sDir = dir([cPath filesep animal filesep 'SpatialDisc' filesep 'Session Data' filesep '*.mat']);
for i = 1:length(sDir)
    load([sDir(i).folder filesep sDir(i).name])
    disp([sDir(i).name ': ' num2str(length(SessionData.Rewarded))]);
end