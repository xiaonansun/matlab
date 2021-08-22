function bhv=behavior_LoadAllSessions(animal)

% This function loads all behavior sessions of an animal.

% animal = 'CSP22';

expTypeStr = 'SpatialDisc';
dataTypeStr = 'Session Data';

if ispc
    rootDir = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
elseif isunix
    rootDir = '/grid/churchland/data/data/Behavior_Simon/';
end

bhvDirCont = dir(fullfile(rootDir,animal,expTypeStr,dataTypeStr));
bhvDirCont=bhvDirCont(~ismember({bhvDirCont.name},{'.','..'})); % Gets rid of the . and .. rows
bhvDirCont=bhvDirCont(find((cell2mat({bhvDirCont.bytes})>1000000))); % Gets rid of small files (few trials)
bhv=struct;
% bhv=cell(1:length(bhvDirCont));

tic
for iFile = 1:length(bhvDirCont)
        bhvFilePath = fullfile(bhvDirCont(iFile).folder,bhvDirCont(iFile).name);
        load(bhvFilePath);
        bhv(iFile).SessionData = SessionData;
        bhv(iFile).SessionName = bhvDirCont(iFile).name(1:end-4);
        disp(['Loaded session: ' bhv(iFile).SessionName])
end
toc
