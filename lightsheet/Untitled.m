exps = twoP_getAcquisitionRecord;
colAnimal = exps{:,1}; colSession = exp{:,6};
if convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('SUNHP')
    baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
elseif convertCharsToStrings(computer) == convertCharsToStrings('GLNXA64') || isunix == 1
    baseDir = '/grid/churchland/data/data/richard_s2p_npy';
end

% subDir = 'logisticRegression';
s2pDir = fullfile('suite2p', 'plane0');
D = cell(size(exps,1),1);
E = struct; E.data = [];
parfor i = 2:size(exps,1)
    try
        animal = coloAnimal{i};
        session = colSession{i};
        %         [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10, false);
        dataPath = fullfile(baseDir,animal,'imaging',session,s2pDir,'data.mat');
        data = load(dataPath ,'data'); data = data.data;
        D{i} = data.neural(1,1,1);
        disp(['Loaded ' animal ' ' session '.']);
    end
end