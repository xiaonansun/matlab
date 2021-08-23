function rateDisc_collectEncodingResults(animal)
% collect results from linear encoding model and collect into single file.
% this is collect explained variance in full, task or motor models.

%% some variables
if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
%     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
%     cPath = 'Q:\BpodImager\Animals\'; %data path on the server
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
dPath = [cPath animal filesep 'SpatialDisc' filesep]; % path for raw data
trainingRange = 'allAudio';

%% block info
load([bPath 'trialInfo_' trainingRange '.mat'], 'recs');
load('allenDorsalMapSM.mat', 'dorsalMaps')
allenMask = dorsalMaps.allenMask;

%% go through recordings and compute variance and prediction error
modelVarMaps = NaN(sum(~allenMask(:)), 3, length(recs), 'single');
for iRecs = 1 : length(recs)
    try
        % load data
        fPath = [dPath filesep recs(iRecs).name filesep]; %session data path
        load([fPath 'mask.mat'], 'mask');
        load([fPath 'opts2.mat'], 'opts');
        load([fPath 'orgfullcorr.mat'], 'fullMap');
        load([fPath 'taskregData.mat'], 'taskMap');
        load([fPath 'motorregData.mat'], 'motorMap');
        
        % align to allen
        fullMap = alignAllenTransIm(arrayShrink(fullMap,mask,'split'),opts.transParams); %align to allen
        taskMap = alignAllenTransIm(arrayShrink(taskMap,mask,'split'),opts.transParams); %align to allen
        motorMap = alignAllenTransIm(arrayShrink(motorMap,mask,'split'),opts.transParams); %align to allen
        
        % larger array
        modelVarMaps(:,1,iRecs) = arrayShrink(fullMap, allenMask, 'merge');
        modelVarMaps(:,2,iRecs) = arrayShrink(taskMap, allenMask, 'merge');
        modelVarMaps(:,3,iRecs) = arrayShrink(motorMap, allenMask, 'merge');
        fprintf('Recording %d/%d finished.\n', iRecs, length(recs));
    end
end

%% save output
save([bPath 'modelMaps_' trainingRange '.mat'],'modelVarMaps');