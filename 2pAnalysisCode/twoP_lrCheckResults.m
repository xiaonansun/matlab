function C=twoP_lrCheckResults

docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';

% Specify session to load
exps = GetGoogleSpreadsheet(docid);

baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
subDir = 'logisticRegression';

C = {};

for i = 2:size(exps,1)
    C{i-1,1} = exps{i,1};
    C{i-1,2} = exps{i,6};
    disp(['Checking session: ' C{i-1,1} ' ' C{i-1,2}]);
    C{i-1,3} = ~isempty(dir([fullfile(baseDir,exps{i,1},'imaging',exps{i,6},subDir) filesep '*_lr.mat']));
    try
        file = dir([fullfile(baseDir,exps{i,1},'imaging',exps{i,6},subDir) filesep '*_lr.mat']);
        load(fullfile(file.folder,file.name),'lr');
        C{i-1,4} = size(lr.iEpoch,2);
    end
    clear lr
end

T = cell2table(C);
writetable(T,[fullfile(baseDir,'logisticRegressionAnalysisSuccess.csv')]);