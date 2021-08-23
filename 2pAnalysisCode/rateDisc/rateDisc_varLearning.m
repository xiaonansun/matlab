function rateDisc_varLearning(fileExt)
% code to compare beta kernels from learning data in different learning 
% episodes and transgenic mice.

if ~exist(fileExt,'var') || isempty(fileExt)
    fileExt = 'org'; %non-orthogonalized model by default
end

% cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
varThresh = 4; %threshold in SDUs for recordings with too high variance (indicative of bad hemo correction)

%% run over animals
% for iAnimals = 1 : length(animals)
for iAnimals = 1 : 8
    
    recCnt = 0;
    cAnimal = animals{iAnimals}; % current animal
    
    % check if animal belongs to a group
    bPath = [cPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    load([bPath 'mask.mat'],'allenMask');
    load([bPath 'trialInfo.mat'], 'recs');
    load([bPath 'sessionVar.mat'],'sessionVar');
    
    useIdx = zscore(sessionVar(end,:)) < varThresh; %only used recordings without excessive variance
    recs(~useIdx) = [];
    recIdx = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs); %identify different training episodes
    fprintf('Rejected %d/%d recordings for excessive variance\n',sum(useIdx),length(useIdx));
    
    for iRecs = 1 : length(recs)
        
        fPath = [cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %session data path
        load([fPath 'Vc.mat'],'bTrials');
        
        if length(bTrials) > 100 && exist([fPath 'dimBeta.mat'],'file') %dont use session with too low trialcount
            
            recCnt = recCnt + 1;
            load([fPath 'opts2.mat'],'opts');
            load([fPath fileExt 'dimBeta.mat'],'dimBeta');
            load([fPath fileExt 'regData.mat'],'regIdx', 'regLabels', 'rejIdx');
            
            load([fPath 'Vc.mat'],'U');
            U = alignAllenTransIm(single(U),opts.transParams); %align to allen
            U = arrayShrink(U, allenMask, 'merge');
            
            if recCnt == 1
                regBeta = cell(1,length(regLabels)); %this will keep a running average over each beta kernel across all recordings
                regMeans = NaN(sum(~allenMask(:)), length(regLabels), length(recs), 'single'); %this will keep a single image for each kernel/recording (average over time)
            end
            
            for iRegs = 1 : length(regLabels)
                if sum(regIdx(~rejIdx) == iRegs) > 0 %regressors are present
                    
                    cData = U * dimBeta(regIdx(~rejIdx) == iRegs, :)';
                    regMeans(:, iRegs, iRecs) = nanmean(cData(:,1:min([sum(regIdx(~rejIdx) == iRegs) ceil(opts.frameRate/2)])),2); %keep average over all regressors
                    cData(isnan(cData(:))) = 0;
                    
                    % compare regressor count in average and current recording
                    if size(regBeta{iRegs},2) < size(cData,2)
                        regBeta{iRegs} = [regBeta{iRegs} zeros(sum(~allenMask(:)), size(cData,2) - size(regBeta{iRegs},2), 'single')]; %add columns for new regressors if needed
                    else
                        cData = [cData zeros(sum(~allenMask(:)),size(regBeta{iRegs},2) - size(cData,2), 'single')]; %add columns for new regressors if neede
                    end
                    regBeta{iRegs} = regBeta{iRegs} + ((cData - regBeta{iRegs})/recCnt); % update mean for current regressors
                    
                end
            end
        end
        fprintf('Recording %d/%d finished.\n', iRecs, length(recs));
        
        %% save beta for current part if new training episode has been reached
        if recIdx(iRecs) < recIdx(iRecs+1) %last recording of current training episode
            save([bPath fileExt 'betaMaps_episode' num2str(recIdx(iRecs)) '.mat'], 'regBeta', 'regMeans', 'recs', 'regLabels', '-v7.3');
            recCnt = 0;
        end
    end
end