function [SessionData,behaviorFilePath] = twoP_loadBehaviorSession(animal,session,bhvFName)
% dbstop 23
minFileSize=1024*8; %minimum size of behavior file (SessionData)

behaviorRootDir= '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
behaviorSubDir = 'SpatialDisc\Session Data';

if length(session) == 6
    sessionDate = datestr(datenum(session,'yymmdd'),'mmmdd_yyyy');
elseif length(session) == 7
    sessionDate = datestr(datenum(session(1:end-1),'yymmdd'),'mmmdd_yyyy');
elseif length(session) == 8
    sessionDate = datestr(datenum(session,'yyyymmdd'),'mmmdd_yyyy');
elseif length(session) == 9
    sessionDate = datestr(datenum(session(1:end-1),'yyyymmdd'),'mmmdd_yyyy');
end

behaviorFileList = dir(fullfile(behaviorRootDir,animal,behaviorSubDir));
dateIdx = find(contains({behaviorFileList.name},sessionDate));
dateIdx(cell2mat({behaviorFileList(dateIdx).bytes}) < minFileSize)=[];

if ~isempty(bhvFName) % load behavior file from a specific file, as specified on the google sheets
    behaviorFilePath = [behaviorFileList(dateIdx(1)).folder filesep bhvFName '.mat'];
    SessionData = load(behaviorFilePath); % loads behavior file that corresponds to the imaging session specified in the google sheets document
else
    if length(session) == 6 || length(session) == 8
        behaviorFilePath = (fullfile(behaviorFileList(dateIdx(1)).folder, behaviorFileList(dateIdx(1)).name));
    elseif (length(session) == 7 || length(session) == 9) && session(end) == 'a'
        behaviorFilePath = (fullfile(behaviorFileList(dateIdx(2)).folder, behaviorFileList(dateIdx(2)).name));
    elseif (length(session) == 7 || length(session) == 9) && session(end) == 'b'
        behaviorFilePath = (fullfile(behaviorFileList(dateIdx(3)).folder, behaviorFileList(dateIdx(3)).name));
    elseif (length(session) == 7 || length(session) == 9) && session(end) == 'c'
        behaviorFilePath = (fullfile(behaviorFileList(dateIdx(4)).folder, behaviorFileList(dateIdx(4)).name));
    end
end

SessionData = load(behaviorFilePath); % loads behavior file that corresponds to the imaging session