function opto_enterNotes(animal)
% Manual checking and entry of notes to indicate location of optogenetics stimulation for a directory of optogenetics
% behavior data

sDir = dir(['Z:\data\Behavior_Simon\' animal '\SpatialDisc\Session Data']);
for i=3:numel(sDir)
    try
    load([sDir(i).folder filesep sDir(i).name]);
    end
    
    error = cell(1,numel(sDir));
    try
        disp(['Session file name: ' sDir(i).name ...
            '; Session file size: ' num2str(sDir(i).bytes/1000000) ...
            'MB; opto duration is: ' num2str(max(SessionData.optoDur)) ...
            '; start time: ' datestr(SessionData.TrialStartTime(1),'yyyymmdd HH:MM:SS')]);
        if isempty(SessionData.Notes{1}) && max(SessionData.optoDur) > 0
            SessionData.Notes{1} = input('Location field is empty, enter location (leave empty if trial should be ignored): ','s');
        else
            disp(['Location is: ' SessionData.Notes{1}])
        end
        save([sDir(1).folder filesep sDir(i).name],'SessionData');
    catch ME
        disp([sDir(i).name ', ' ME.identifier ': ' ME.message]);
        error{i} = ME;
    end
end
