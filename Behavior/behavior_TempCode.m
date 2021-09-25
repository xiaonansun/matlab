
exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1); colSession = exps(:,6);

S = twoP_settings;
bhvPlotDir = fullfile(S.dir.imagingRootDir,'behavior_plots');
if ~exist(bhvPlotDir,'dir')
    mkdir(bhvPlotDir);
end

for i = 2:size(exps,1)
    try
    animal = colAnimal{i};
    session = colSession{i};
%     [~,bhv]=twoP_loadImgBhvData(animal,session, true, 10, true);
    bhvDir = fullfile(S.dir.bhvRootDir,animal,S.dir.bhvSubDir);
    [bhv,~] = twoP_loadBehaviorSession(animal,session,bhvDir);
    figA = figure(1);
    RateDisc_getDiscPerformance(bhv)
    exportgraphics(gca,fullfile(bhvPlotDir,[animal '_' session '.pdf']))
    close all
    end
end