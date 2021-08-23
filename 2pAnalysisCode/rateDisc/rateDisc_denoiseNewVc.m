function rateDisc_denoiseNewVc(animals)

if ispc
    cPath = 'Y:\data\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals';
end

% dataOverview = rateDiscRecordings;
opts.nSVD = 200;
opts.frameRate = 30;
opts.blockDims = 500; %number of dimensions from SVD per block
opts.maxLag = 5; %lag for autocorrelation
opts.autoConfidence = 0.99; %confidence for autocorrelation test
opts.autoThresh = 1.5; %threshold for autocorrelation test
opts.snrThresh = 1.6; %threshold for SNR test

% animals = unique(dataOverview(:,1));
animals = unique(animals);
for iAnimal = 1 : length(animals)
    
    recs = ls([cPath animals{1} filesep 'SpatialDisc' filesep]);
    recs = recs(~ismember(recs(:,1), '.'), :);
    
    for iRecs = 1 : size(recs,1)
        try
            fPath = [cPath animals{1} filesep 'SpatialDisc' filesep strtrim(recs(iRecs,:)) filesep];
            load([fPath filesep 'Vc.mat'], 'Vc');
            
            if size(Vc,1) > 200
                load([fPath filesep 'blueV.mat']);
                load([fPath filesep 'hemoV'], 'hemoV');
                
                U = U(:, :, 1:200);
                blueV = blueV(1:200, :, :);
                hemoV = hemoV(1:200, :, :);
                
                [Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
                save([fPath 'Vc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
                save([fPath 'HemoCorrection.mat'],'regC','T', 'hemoVar')
            end
            fprintf('Done: %s\n', fPath)
        
        catch
            fprintf('!! Warning: Could not re-compute Vc in %s !!\n', fPath)
        end
    end
end

