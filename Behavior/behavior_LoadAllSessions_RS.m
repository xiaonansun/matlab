function bhv=behavior_LoadAllSessions(animal,bhvRootDir,bhvFilePathList)

% This function loads all behavior sessions of an animal.
% Inputs
% 1. animal: can be a single char or string variable or a cell with
% multiple animals
% 2. bhvRootDir: user manual of the root directory containing behavior files. The
% subdirectory of the root should be animal IDs
% 3. bhvFilePathList (optional): struct list of file paths generated from dir()
% function

% animal = 'CSP22';

expTypeStr = 'SpatialDisc';
dataTypeStr = 'Session Data';

if ispc
    rootDir = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
elseif isunix
    rootDir = '/grid/churchland/data/data/Behavior_Simon/';
elseif ismac
    if ~exist('bhvRooDir','var') || isempty(bhvRootDir)
    rootDir = '/Users/xiaonansun/Documents/data/twoP';
    subDir = 'SpatialDisc';
    end
end

% if ispc || isunix
% bhvDirCont = dir(fullfile(rootDir,animal,expTypeStr,dataTypeStr));
% bhvDirCont=bhvDirCont(~ismember({bhvDirCont.name},{'.','..'})); % Gets rid of the . and .. rows
% bhvDirCont=bhvDirCont(find((cell2mat({bhvDirCont.bytes})>1000000))); % Gets rid of small files (few trials)
% bhv=struct;
% end

if exist('bhvFilePathList','var') || ~isempty(bhvFilePathList)
%     bhvDirCont=bhvDirCont(~ismember({bhvDirCont.name},{'.','..'})); % Gets rid of the . and .. rows
    bhvDirCont=bhvFilePathList;
%     bhvDirCont=bhvDirCont(find((cell2mat({bhvDirCont.bytes})>500000))); % Gets rid of small files (few trials)
    bhv=struct;
end

% bhv=cell(1:length(bhvDirCont));
%% Combine behavior control data from individual sessions into a single struct

tic
for iFile = 1:length(bhvDirCont)
        bhvFilePath = fullfile(bhvDirCont(iFile).folder,bhvDirCont(iFile).name);
        load(bhvFilePath);
%         bhv(iFile).SessionData = SessionData;
        [file_path,~,~] = fileparts(fullfile(bhvDirCont(iFile).folder,bhvDirCont(iFile).name));
        idxFS = regexp(file_path,filesep);
        session = file_path(idxFS(end)+1:end);
        bhv(iFile).SessionData = cBhv;
        bhv(iFile).SessionName = session;
        disp(['Loaded animal: ' animal ' and session: ' bhv(iFile).SessionName])
end
toc
