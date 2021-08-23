function rateDisc_areaCorr(group,trainingRange)

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

savePath = 'R:\BpodImager\areaCorr\'; %local save path
load('allenDorsalMapSM.mat', 'dorsalMaps')
allenMask = dorsalMaps.allenMask;
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
animals = animals(contains(animals,group));

%get area mask
targArea = 'RSPd'; %retrosplinal
areaIdx = dorsalMaps.edgeOutline{strcmpi(dorsalMaps.labels, targArea)};
areaMask = poly2mask(areaIdx(:,2), areaIdx(:,1), size(allenMask,1), size(allenMask,2));
areaMask = arrayShrink(areaMask, allenMask);

%% load data and do some analysis
figure('Renderer','painters');
allCorr = NaN(sum(~allenMask(:)), length(animals), 'single');
for iAnimals = 1:length(animals)
    
    recs = rateDisc_getRecsForAnimal(animals{iAnimals},trainingRange);% end
    areaCorr = NaN(sum(~allenMask(:)), length(recs), 'single');
    for iRecs = 1 : length(recs)
        try
            fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep];
            try
                load([fPath 'rsVc.mat'], 'Vc');
            catch
                load([fPath 'Vc.mat'], 'Vc');
            end
            load([fPath 'alignU.mat'], 'U');
            
            U = rateDisc_removeOutline(U, 10);
            U = arrayShrink(U, allenMask, 'merge');
            Vc = reshape(Vc, size(Vc,1), []);
            Vc(:,isnan(Vc(1,:))) = [];
            
            covV = cov(Vc');
            varP = dot((U*covV)', U'); % 1 x P
            
            covP = mean(U(areaMask,:)*covV,1)*U'; % 1 x P
            stdPxPy = mean(varP(areaMask)).^0.5 * varP.^0.5; % 1 x P
            areaCorr(:,iRecs) = covP./stdPxPy; % 1 x P
            fprintf('Done: %s - %s !!\n', animals{iAnimals}, recs(iRecs).name)
            
        catch ME
            fprintf('!! Warning: Failed to run model %s - %s !!\n', animals{iAnimals}, recs(iRecs).name)
            disp(ME.message);
        end
    end
    
    %% save result
    if ~exist(savePath,'dir')
        mkdir(savePath);
    end
    save([savePath 'areaCorr_' animals{iAnimals} '_' targArea], 'areaCorr');
    allCorr(:,iAnimals) = nanmean(areaCorr,2);
    
    %% show figure
    subplot(1,length(animals),iAnimals);
    imageScale(arrayShrink(allCorr(:,iAnimals), allenMask, 'split'));
    colormap(colormap_blueblackred(256)); caxis([-1 1]);
    title(animals{iAnimals}); drawnow;
end

%% show figure
figure('Renderer','painters');
imageScale(arrayShrink(nanmean(allCorr,2).^2, allenMask, 'split'));
colormap(colormap_blueblackred(256)); caxis([-1 1]); drawnow;
title([group ' - all']);

