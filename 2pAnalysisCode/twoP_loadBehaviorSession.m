function [bhv,behaviorFilePath] = twoP_loadBehaviorSession(animal,session,bhvDir)
% dbstop 23

% --- Load the google sheets document "2photon acquisition record" --- %
docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';
expTable=GetGoogleSpreadsheet(docid); % this function (GetGoogleSpreadsheet.m) needs to be downloaded
bhvColIdx=find(contains(expTable(1,:),'Behavior file name'));
iFolderColIdx=find(contains(expTable(1,:),'Folder'));

try
bhvRowIdx = find(contains(expTable(:,iFolderColIdx),session));
bhvFName = expTable{bhvRowIdx(contains(expTable(bhvRowIdx),animal)),bhvColIdx};
catch ME
    disp([ME.message]);
    disp('Cannot load session. Please check the session name input.');
    analysis.error.behaviorTable = ME.message;
end

minFileSize=1024*8; %minimum size of behavior file (SessionData)

% behaviorRootDir= '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
% behaviorSubDir = 'SpatialDisc\Session Data';

if length(session) == 6
    sessionDate = datestr(datenum(session,'yymmdd'),'mmmdd_yyyy');
elseif length(session) == 7
    sessionDate = datestr(datenum(session(1:end-1),'yymmdd'),'mmmdd_yyyy');
elseif length(session) == 8
    sessionDate = datestr(datenum(session,'yyyymmdd'),'mmmdd_yyyy');
elseif length(session) == 9
    sessionDate = datestr(datenum(session(1:end-1),'yyyymmdd'),'mmmdd_yyyy');
end

behaviorFileList = dir(bhvDir);
dateIdx = find(contains({behaviorFileList.name},sessionDate));
dateIdx(cell2mat({behaviorFileList(dateIdx).bytes}) < minFileSize)=[];

if ~isempty(bhvFName) % load behavior file from a specific file, as specified on the google sheets
    behaviorFilePath = fullfile(behaviorFileList(dateIdx(1)).folder, bhvFName '.mat');

    disp('The corresponding Bpod behavior data has been loaded.')
    return
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

load(behaviorFilePath,'SessionData'); % loads behavior file that corresponds to the imaging session specified in the google sheets document
bhv = SessionData;