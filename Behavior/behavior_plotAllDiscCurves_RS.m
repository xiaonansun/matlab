function behavior_plotAllDiscCurves(animal,bhvFilePathList)
% generate behavior plot
% 2023-06-27 edited by Richard Sun

%%
if ~exist('bhvFileList','var') && isempty(bhvFilePathList)
%     animal = 'Fez72'; % uncomment this line to run original code
    sDir = dir(['Z:\data\Behavior_Simon\' animal '\SpatialDisc\Session Data']);
else
    sDir = bhvFilePathList;
end

%%
% for i=3:numel(sDir)
for i=1:numel(sDir)
    %%
    if sDir(i).bytes > 1000
    try
    load([sDir(i).folder filesep sDir(i).name]);
    SessionData = cBhv; % RS edit
    catch ME
    end
    error = cell(1,numel(sDir));
    try
        [file_path,~,~] = fileparts(fullfile(sDir(i).folder,sDir(i).name));
        idxFS = regexp(file_path,filesep);
        session = file_path(idxFS(end)+1:end);
        disp(['Session file name: ' session ...
            '; Session file size: ' num2str(sDir(i).bytes/1000000) ...
            'MB; start time: ' datestr(SessionData.TrialStartTime(1),'yyyymmdd HH:MM:SS')]);
        cInd=1:SessionData.nTrials;
        [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(SessionData, cInd);
        hDisc=figure(1);
        plot(distRatio,pChoseHigh,'-k');
        title([animal ' ' session]); xlabel('Fraction of right-sided events'); ylabel('Fraction of right-sided choice');
        fig_configAxis(gca);
        pause;
        close(hDisc); clear SessionData;
    catch ME
        disp([sDir(i).name ', ' ME.identifier ': ' ME.message]);
        error{i} = ME;
    end
    end
end

