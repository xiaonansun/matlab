cPath = 'Y:\data\BpodImager\Animals\Plex02\SpatialDisc\';

recs = dir(cPath);
allTrials = cell(1, length(recs));
for iRecs = 1 : length(recs)
    
    try 
        load([cPath recs(iRecs).name '\Vc.mat'],'bTrials')
        allTrials{iRecs} = bTrials;
        
    catch
        allTrials{iRecs} = [];
    end
end

%% do plot
figure; hold on
for iRecs = 1 : length(recs)
    if ~isempty(allTrials{iRecs})
        plot(allTrials{iRecs})
    else disp(iRecs);
    end
end