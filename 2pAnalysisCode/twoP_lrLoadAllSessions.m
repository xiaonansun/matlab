function LR=twoP_lrLoadAllSessions(loadCombined)
% lr is single session logistic regression data
% LR is a struct that contains the lr struct from all sessions
% loadCombined = 0 looks through all all session folders and creates a new
% LR struct
% loadCombined = 1 looks for and loads the LR struct that has been save previously

if ~exist('loadCombined','var') || isempty(loadCombined)
    loadCombined = 1;
end

% Specify session to load
exps = twoP_getAcquisitionRecord;

baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
subDir = 'logisticRegression';
s2pDir = 'suite2p\plane0';
LR = struct;
iStart = 2;

% animalList = {'CSP27'; 'CSP29'; 'CSP30'; 'Plex50'; 'Plex51'; 'Fez51'; 'Fez57'; 'Fez59'};
if loadCombined == 1
    D = dir(fullfile(baseDir,'*logisticRegressionAllAnimals.mat'));
    DD = zeros(length(D),1);
    for i = 1:length(D); DD(i) = datenum(D(i).date); end
    [~,j]=max(DD);
    disp('Loading previously combined logistic regression analysis data...');
    load(fullfile(D(j).folder,D(j).name),'LR');
    disp(['Loaded ' fullfile(D(j).folder,D(j).name)]);
    disp('Done!');
elseif loadCombined == 0
    disp('Loading and combining logistic regression analysis data from individual sessions...');
    for i = iStart:size(exps,1)
        LR(i).animal = exps{i,1};
        LR(i).session = exps{i,6};
        LR(i).date = datetime(exps{i,2},'InputFormat','MM/dd/yyyy');
        LR(i).area = exps{i,3};
        LR(i).depth = exps{i,4};
        LR(i).cvAcc =  []; LR(i).cvAcc_r = []; LR(i).mcvAcc_nr_rep = [];
        LR(i).iEpoch = [];
        LR(i).bMaps = [];
        LR(i).allAUC = [];
        LR(i).segFrames = [];
        LR(i).shufMu = [];
        LR(i).shufSigma = [];
        try
            load([fullfile(baseDir,LR(i).animal,'imaging',LR(i).session,subDir) filesep LR(i).animal '_' LR(i).session '_lr.mat'],'lr');
            load([fullfile(baseDir,LR(i).animal,'imaging',LR(i).session,s2pDir,'data.mat')],'data');
            %     load([fullfile(baseDir,animal,'imaging',session) filesep s2pDir filesep 'data.mat']);
            LR(i).cvAcc = lr.cvAcc; LR(i).cvAcc_r = lr.cvAcc_r; LR(i).mcvAcc_nr_rep = lr.mcvAcc_nr_rep;
            LR(i).bMaps = lr.bMaps; LR(i).allAUC = lr.allAUC; LR(i).iEpoch = lr.iEpoch; LR(i).segFrames = lr.segFrames;
            LR(i).Vc = lr.Vc;
            LR(i).trialNumbers = data.trialNumbers;
            LR(i).idx_redcell = data.idx_redcell; LR(i).idx_notredcell = data.idx_notredcell;
            mu = zeros(size(lr.shufAUC,1),size(lr.shufAUC,2)); sigma = zeros(size(lr.shufAUC,1),size(lr.shufAUC,2));
            shufAUC = lr.shufAUC;
            for x = 1:size(shufAUC,1)
                parfor y = 1:size(shufAUC,2)
                    [mu(x,y),sigma(x,y)] = normfit(squeeze(shufAUC(x,y,:)));
                end
            end
            LR(i).shufMu = mu; LR(i).shufSigma = sigma;
        end
        disp([LR(i).animal ' ' LR(i).session ' loaded.']);
    end
    disp('Completed: loading and combining logistic regression analysis data.');
    save(fullfile(baseDir,[datestr(now,'yyyy-mm-dd_HHMMSS') '_logisticRegressionAllAnimals.mat']),'LR');
    disp(['Saved combined logistic regression analysis data as: ' fullfile(baseDir,[datestr(now,'yyyy-mm-dd_HHMMSS') '_logisticRegressionAllAnimals.mat'])]);
end