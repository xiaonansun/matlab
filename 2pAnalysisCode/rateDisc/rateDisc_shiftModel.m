function rateDisc_shiftModel(animal,cRec,cShift)

if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %data path on the server
%     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
%     cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
end
maxShift = 45; %determines the range in which data will be shifted. abs(cShift) cannot be > maxShift.
ridgeFolds = 10;
if abs(cShift) > maxShift; error(['cShift cannot be larger than maxShift (' num2str(maxShift) ')']); end

% try
    fPath = [cPath animal filesep 'SpatialDisc' filesep cRec filesep];
    load([fPath 'Vc.mat'], 'U');
    load([fPath 'mask.mat'], 'mask');
    load([fPath 'opts2.mat'], 'opts');
    
    load([fPath 'interpVc.mat'], 'Vc', 'frames');
    load([fPath 'regData.mat'], 'fullR');
    
    load('allenDorsalMapSM.mat', 'dorsalMaps')
    allenMask = dorsalMaps.allenMask;
    
    % crop and shift Vc
    Vc = reshape(Vc, size(Vc,1), frames, []);
    Vc(:,1:maxShift,:) = 0;
    Vc(:,end - maxShift + 1 : end,:) = 0;
    Vc = circshift(Vc,cShift,2);
    Vc = reshape(Vc, size(Vc,1), []);
    
    % run cross-validation
    rng default %reset randum number generator
    randIdx = randperm(size(Vc,2)); %generate randum number index
    foldCnt = floor(size(Vc,2) / ridgeFolds);
    Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
    
    for iFolds = 1:ridgeFolds
        
        dataIdx = true(1,size(Vc,2));
        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        if iFolds == 1
            [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', fullR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
        else
            [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', fullR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
        end
        Vm(:,~dataIdx) = (fullR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data
        
        if rem(iFolds,ridgeFolds/5) == 0
            fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
        end
    end
    
    corrMat = rateDisc_modelCorr(Vc,Vm,U) .^2; %compute explained variance after shifting data by one frame
    corrMat = arrayShrink(corrMat,mask,'split');
    corrMat = alignAllenTransIm(corrMat,opts.transParams);
    corrMat = nanmean(arrayShrink(corrMat(:,1:size(allenMask,2)),allenMask, 'merge'));
    
    tPath = [fPath 'shiftVals' filesep];
    if ~exist(tPath, 'dir')
        mkdir(tPath);
    end
    
    save([tPath 'shiftCorr_' num2str(cShift) '.mat'], 'corrMat', 'maxShift', 'ridgeFolds');
    fprintf('Done: %s - %s\n', animal, cRec)
% catch
%     fprintf('!! Warning: Failed to run model %s !!\n', recs(iRecs).name)
% end