function twoP_lrOptimizeDecoder

docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';

% Specify session to load
exps = GetGoogleSpreadsheet(docid);
iStart = 2;

for i = iStart:size(exps,1)
    close all
    
    animal = exps{i,1};
    session = exps{i,6};
    
    disp(['Processing animal ' animal ' and session ' session '.']);
    
    try
        baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy'; % directory of all 2p data
        lrDir = fullfile(baseDir,animal,'imaging',session,'logisticRegression'); % directory to save the current logistic regression analysis data and figures
        
        if exist(fullfile(lrDir,[animal '_' session '_LogisticRegressionOptimization.mat']),'file')
            disp('Optimization analyses have been completed for this session.');
            continue
        end

        [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10);
        
        % ----- Makes adjustments to the twoP data struct ----- %%
        data = twoP_adjustData(data,SessionData);
        
        
        if ~exist(lrDir,'dir') % if the current analysis directory doesn't exist, then create one
            mkdir(lrDir);
        end
        lrDataPath = fullfile(lrDir, [animal '_' session '_lr.mat']); % create a file path for the current session's data (.mat file)
        
        % Initialize some parameters for logistic regression
        regType = 'lasso'; %lasso or ridge
        stepSize = [];
        decType = 'allChoice';
        segIdx = [1 0.75 1.25 0.5 1];
        reps = 10; % number of times to randomly subsample non-red cells
        minReps = 10; % number of minimum subsamples for data to be included for analysis
        opts.preStim = data.trialStimFrame*data.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
        opts.frameRate = 1000/data.msPerFrame; % Frame rate of imaging
        sRate = opts.frameRate;
        segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment
        loadExistingData = false;
        cBhv = selectBehaviorTrials(SessionData, data.trialNumbers); %% very important in matching trial indices
        
        if exist(lrDataPath,'file') && loadExistingData == true
            disp('Existing analyzed data exists, loading...');
            load(lrDataPath)
            if size(lr.cvAcc_nr_rep,1) >= minReps
                subSampling = false;
            else
                disp('There are less than 10 subsamples, will re-run...');
                subSampling = true;
            end
        end
        
        if ~exist('subSampling','var') || isempty(decType)
            subSampling = false;
        end
        
        % Aligns imaging data based on the behavior data
        Vc = data.neural;
        Vc_r = data.neural(data.idx_redcell ,:,:);
        Vc_nr = data.neural(data.idx_notredcell,:,:);
        
        % Define the useTrials parameter
        % useTrials = floor(size(Vc,3)*0.8);
        % if length(data.trialNumbers) >= 500
        %     useTrials = 400;
        % elseif length(data.trialNumbers) >= 400 && length(data.trialNumbers) < 500
        %     useTrials = 300;
        % elseif length(data.trialNumbers) >= 300 && length(data.trialNumbers) < 400
        %     useTrials = 200;
        % elseif length(data.trialNumbers) < 200
        %     disp('There are fewer than 200 trials in this session, terminating analysis...');
        %     return
        % end
        
        Vc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts); %aligned to different trial episodes
        Vc_r = rateDisc_getBhvRealignment(Vc_r, cBhv, segFrames, opts); %aligned to different trial episodes
        Vc_nr = rateDisc_getBhvRealignment(Vc_nr, cBhv, segFrames, opts); %aligned to different trial episodes
        
        cvAcc = nan(size(Vc,3),size(Vc,2), 'single');
        parfor i = 1:size(Vc,3)
            [cvAcc(i,:), bMaps, betaNeuron, mdlAll, trialCnt, allAUC] = rateDisc_logDecoder(Vc, [], cBhv, i, 0, regType, stepSize, decType);
            disp(['Current trial ' num2str(i)]);
        end
        
        imagesc(cvAcc); colormap(parula); cb = colorbar;
        ylabel('useTrials');
        xlabel('Frame');
        title([animal ' ' session]);
        save(fullfile(lrDir,[animal '_' session '_LogisticRegressionOptimization.mat']),'cvAcc');
        exportgraphics(gcf,fullfile(lrDir,[animal '_' session '_LogisticRegressionOptimization.pdf']));
    end
end