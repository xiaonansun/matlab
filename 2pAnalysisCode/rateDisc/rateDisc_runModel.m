function rateDisc_runModel(animal,idx,trainingRange)
%code to run the regression model over a specific recording for an animal
%in the 'rateDiscRecordings' file. 'iRecs' is the index of the recording
%out of the 'trainingRange' selection. If trainingRange is undefined will
%default to 'audioDisc' for auditory discrimination sessions.

if ispc
    %     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
    %     cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
    cPath = 'Q:\BpodImager\Animals\'; %local data path
else
    %     cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/'; %path to grid
end

if ~exist('trainingRange','var') || isempty(trainingRange)
    trainingRange = 'audioDisc'; %use this to run analysis only in a certain range of training
end

%% get list of recordings and reference image
% try
%     load([cPath animal filesep 'blockData' filesep 'trialInfo.mat'],'recs')
% catch
%     disp('Couldnt find trialInfo in blockData folder. Moving all recs instead.');
recs = rateDisc_getRecsForAnimal(animal,trainingRange);% end

% run over all recordings in selected range
if idx == 0
    idx = 1 : length(recs);
end

for iRecs = idx
    try
    %     fPath = [cPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep];
    %     load([fPath 'Vc.mat'], 'U');
    %     load([fPath 'opts.mat'], 'opts');
    
    %     if size(Vc,1) > 200
    %         load([fPath filesep 'blueV.mat'],'blueV','blueFrametimes','Sv','totalVar','trials','bTrials');
    %         load([fPath filesep 'hemoV'], 'hemoV');
    %
    %         U = U(:, :, 1:200);
    %         blueV = blueV(1:200, :, :);
    %         hemoV = hemoV(1:200, :, :);
    %
    %         [Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
    %         save([fPath 'Vc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
    %         save([fPath 'HemoCorrection.mat'],'regC','T', 'hemoVar')
    %     end
    %     rateDisc_RegressModel(cPath, animal, recs(iRecs).name, 'Widefield');
    %     rateDisc_choiceModel(cPath, animal, recs(iRecs).name, 'Widefield');
    %     rateDisc_choiceModelDLC(cPath, animal, recs(iRecs).name, 'Widefield');
    rateDisc_logRegress(cPath, cPath, animal, recs(iRecs).name, 'lasso', 0, 300);
    fprintf('Done: %s - %d - %s\n', animal, iRecs, recs(iRecs).name)

    catch ME
        fprintf('!! Warning: Failed to run model %s - %s !!\n', animal, recs(iRecs).name)
        disp(ME.message);
    end
end