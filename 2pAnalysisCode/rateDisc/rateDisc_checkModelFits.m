rateDisc_checkModelFits

[dataOverview, ~, ~, ~, segIdx, segLabels, segIdxRealign, cPath, trainDates] = rateDiscRecordings;
targMod = 0; %modality (0 = all trials, 2 = audio, 4 = tactile, 5 = audiotactile)
regType = 'lasso'; %lasso or ridge
tPath = 'X:\smusall\BpodImager\Animals\';
trainDur = 'all'; %period from which imaging data should be used
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask; %mask that is used for all datasets
segFrames = floor(segIdx * opts.frameRate); %max nr of frames per segment

for iAnimals = 1 : 4

    recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
    %get training range of interest
    if ~strcmpi(trainDur, 'all')
        cDates = datenum(cTrainDates{iAnimals}(contains(cTrainLabels, trainDur)));
    else
        cDates = [-inf; inf];
    end
    for iRecs = 1 : length(recs)
        try
            useIdx(iRecs) = datenum(recs(iRecs).name) >= min(cDates) && datenum(recs(iRecs).name) <= max(cDates);
        catch
            useIdx(iRecs) = false;
        end
    end
    recs = recs(useIdx);
    
    
    
    
    Cnt = 1;
    for iRecs = 3 : size(recs,1)
        try
            fPath = [cPath animal filesep 'SpatialDisc' filesep rec filesep]; %Widefield data path

            
            
            fprintf('Done: %s - %s\n', animal, recs(iRecs).name)
        catch
            fprintf('!! Warning: Failed to run model %s !!\n', recs(iRecs).name)
        end
    end
end