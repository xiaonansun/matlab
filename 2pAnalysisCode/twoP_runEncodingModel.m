function twoP_runEncodingModel(animal,session,dType,cPath,Rec)
%% Use this script to run the encoding model on a single session.
animal = 'CSP27'; session = '20200308'; dType = 'twoP'; cPath = 'H:\twoP'; Rec = [];

S = twoP_settings;
[~,behaviorFilePath] = twoP_loadBehaviorSession(animal,session);
twoPdataDir= fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir);
twoPdataPath = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'data.mat');
bhvFilename = regexp(behaviorFilePath,filesep,'split'); bhvFilename = bhvFilename{end};
bhvFilename = regexp(bhvFilename,'\.','split'); bhvFilename=bhvFilename{1}; % need to change the backslash "\" to filesep
bhvVidDir = fullfile(S.dir.bhvVidRootDir,animal,S.dir.bhvVidSubDir,bhvFilename);
if ~exist(bhvVidDir,'dir'); disp('Behavior video folder does not exist, is the behavior data split?'); end

myDirs.bhvFilePath = behaviorFilePath;
myDirs.twoPdataPath = twoPdataPath;
myDirs.twoPdataDir= twoPdataDir;
myDirs.bhvVidDir = bhvVidDir;

%
rateDisc_choiceModel(cPath,animal,session,dType,myDirs);