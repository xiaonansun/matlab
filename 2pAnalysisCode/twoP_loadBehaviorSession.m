function [bhv,behaviorFilePath] = twoP_loadBehaviorSession(animal,session,varargin)
% Input variables
% (1) animal
% (2) session
% (3) data (optional)

if nargin > 2
    data = varargin{3};
end
bhvFName = [];
% dbstop 23

% --- Load the google sheets document "2photon acquisition record" --- %
% docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';
% expTable=GetGoogleSpreadsheet(docid); % this function (GetGoogleSpreadsheet.m) needs to be downloaded
S = twoP_settings;
bhvDir = fullfile(S.dir.bhvRootDir,animal,S.dir.bhvSubDir);
expTable = S.exps;
bhvColIdx=contains(expTable(1,:),'Behavior file name');
iFolderColIdx=contains(expTable(1,:),'Folder');

% idxRow = twoP_findRowAcquisitionRecord(animal,session);
% sFilenames = regexp(expTable{idxRow,9},'\;','split','once');

try
bhvRowIdx = find(contains(expTable(:,iFolderColIdx),session));
bhvFName = expTable{bhvRowIdx(contains(expTable(bhvRowIdx),animal)),bhvColIdx};
catch ME
    disp([ME.message]);
    disp('Behavior file not specified in 2P acquisition record, will look for file in the behavior directory.');
%     analysis.error.behaviorTable = ME.message;
end

minFileSize=1024*8; %minimum size of behavior file (SessionData)

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
    sFilenames = regexp(bhvFName,'\;','split','once');
    if length(sFilenames) > 1
        if ~exist('data','var') || isempty(data)
            load(fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'data.mat'),'data');
        end
        [SessionData,~] = twoP_RepairSession(animal,sFilenames,data.trialCodes);
        behaviorFilePath = fullfile(bhvDir,[strcat(sFilenames{:}) '.mat']);
        save(behaviorFilePath,'SessionData');
        disp(['Behavior data repaired/concatenated, saved as: ' behaviorFilePath]);
    else
        behaviorFilePath = fullfile(behaviorFileList(dateIdx(1)).folder, [bhvFName '.mat']);
    end
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

if ~iscell(behaviorFilePath)
    SessionData = load(behaviorFilePath,'SessionData'); % loads behavior file that corresponds to the imaging session specified in the google sheets document
    bhv = SessionData.SessionData;
end