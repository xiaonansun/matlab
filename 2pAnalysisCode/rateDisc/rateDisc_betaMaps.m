function rateDisc_betaMaps(fileExt)
% code to compare beta kernels from learning data in different learning
% episodes and transgenic mice.

if ~exist('fileExt','var') || isempty(fileExt)
    fileExt = 'org'; %non-orthogonalized model by default
%     fileExt = 'choice'; %non-orthogonalized model by default
end

% cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
% cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
cPath = 'Q:\BpodImager\Animals\'; %path to churchlandNAS
tPath = 'Q:\BpodImager\Animals\'; %path to churchlandNAS
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
trainingRange = 'audioDisc'; %use this to run analysis only in a certain range of training

%% run over animals
% for iAnimals = 1
for iAnimals = 1 : length(animals)
    
    cAnimal = animals{iAnimals}; % current animal
    bPath = [cPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    tbPath = [tPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    if ~exist(tbPath,'dir')
        mkdir(tbPath);
    end
    load([bPath 'mask_allAudio.mat'],'allenMask');
    recs = rateDisc_getRecsForAnimal(cAnimal, trainingRange, cPath);
    load([cPath cAnimal filesep 'SpatialDisc' filesep recs(end).name filesep fileExt 'regData.mat'], 'regLabels'); %get regLabels
    
    clear allU
    fprintf('Next animal: %s\n', cAnimal);
    for iRegs = 1 : length(regLabels)
        
        firstRec = true; %flag to create new 'cBeta'
        for iRecs = 1 : length(recs)
            
            fPath = [cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %session data path
            %there can be NaNs in recs. this is to keep the code from crashing.
            bTrials = []; 
            if exist([fPath 'rsVc.mat'],'file')
                vcFile = 'rsVc.mat';
            else
                vcFile = 'Vc.mat';
            end
            load([fPath vcFile],'bTrials');

            if length(bTrials) > 100 && exist([fPath fileExt 'regData.mat'],'file') %dont use session with too low trialcount
                
                load([fPath 'opts2.mat'],'opts');
                load([fPath fileExt 'dimBeta.mat'],'dimBeta');
                load([fPath fileExt 'regData.mat'],'regIdx', 'regLabels', 'rejIdx');
                cIdx = regIdx == iRegs; %current regressors
                
                if firstRec
                    cBeta = NaN(sum(~allenMask(:)), sum(cIdx), length(recs), 'single'); %initialize cBeta for this regressor
                    firstRec = false;
                end
            
                % get U for first regressor and save to use it for other regressors
                if iRegs == 1
                    load([fPath vcFile],'U');
                    U = alignAllenTransIm(single(U),opts.transParams); %align to allen
                    U = rateDisc_removeOutline(U, 10); % slightly crop outline to avoid rotation artifacts
                    allU{iRecs} = arrayShrink(U, allenMask, 'merge'); clear U
                end
                
                cData = allU{iRecs} * dimBeta(regIdx(~rejIdx) == iRegs, :)';
                if ~isempty(cData)
                    cData = arrayShrink(cData,allenMask,'split');
                    cBeta(:, ~rejIdx(cIdx), iRecs) = arrayShrink(cData, allenMask, 'merge'); 
                    clear cData
                end
            else
                if iRegs == 1 %give feedback on first regressor
                    fprintf('Ignored recording %d for insufficient trials or no modeling results\n', iRecs);
                end
            end
        end
        cBeta = nanmean(cBeta,3); %average over all recordings
        save([tbPath fileExt 'betaMaps_' trainingRange '_' regLabels{iRegs} '.mat'], 'cBeta', 'recs', 'regLabels', '-v7.3');
        fprintf('Finished regressor: %s for animal %s\n', regLabels{iRegs}, cAnimal);
    end
    fprintf('%s - All done!\n', cAnimal);
end
