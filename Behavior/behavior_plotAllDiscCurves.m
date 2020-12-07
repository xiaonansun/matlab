function behavior_plotAllDiscCurves(animal)
% behavior plot
animal = 'Fez72';

sDir = dir(['Z:\data\Behavior_Simon\' animal '\SpatialDisc\Session Data']);
for i=3:numel(sDir)
    if sDir(i).bytes > 1000
    try
    load([sDir(i).folder filesep sDir(i).name]);
    catch ME
    end
    error = cell(1,numel(sDir));
    try
        disp(['Session file name: ' sDir(i).name ...
            '; Session file size: ' num2str(sDir(i).bytes/1000000) ...
            'MB; start time: ' datestr(SessionData.TrialStartTime(1),'yyyymmdd HH:MM:SS')]);
        cInd=1:SessionData.nTrials;
        [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(SessionData, cInd);
        hDisc=figure(1);
        plot(distRatio,pChoseHigh,'-k');
%         title(['']);
        pause;
        close(hDisc); clear SessionData;
    catch ME
        disp([sDir(i).name ', ' ME.identifier ': ' ME.message]);
        error{i} = ME;
    end
    end
end

