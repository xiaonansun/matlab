function rateDisc_evaluateSessionSVD(animal)
%% some variables
if ispc
%     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
dPath = [cPath animal filesep 'SpatialDisc' filesep]; % path for raw data
nDims = 200; %maximal nr of local dimensions

load([bPath 'trialInfo_allAudio.mat'],'recs');
load([bPath 'mask_allAudio.mat'],'allenMask','mask','xRange','yRange');
mask = allenMask(yRange,xRange) | mask;
[~, ~, ~, ~, segIdx] = rateDiscRecordings;
trainingRange = 'allAudio';

%% go through recordings and compute variance and prediction error
tic;
sessionVar = NaN(nDims, length(recs));
sessionVarMap = NaN(sum(~allenMask(:)), length(recs), 'single');
sessionAvgVarMap = NaN(sum(~allenMask(:)), length(recs), 'single');
for iRecs = 1 : length(recs)

    fPath = [dPath recs(iRecs).name filesep]; disp(fPath); %session data path
    load([fPath 'opts2.mat'],'opts');
    if exist([fPath 'rsVc.mat'], 'file')
        load([fPath 'rsVc.mat'],'Vc', 'U', 'bTrials'); %use downsampled Vc
        opts.frameRate = 15; %adjust framerate
    else
        load([fPath 'Vc.mat'],'Vc', 'U', 'bTrials');
    end
    U = removeEdges(U, 5, NaN); %make sure there are no edge effects
    U = alignAllenTransIm(single(U),opts.transParams); %align to allen
    
    % get behavioral data
    bhvFile = dir([fPath animal '*SpatialDisc*.mat']);
    load([fPath bhvFile.name], 'SessionData');
    bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dat
    
    % make sure Vc is zero-mean
    [~, B, C] = size(Vc);
    Vc = reshape(Vc(1:nDims,:,:), nDims,[]);
    Vc = bsxfun(@minus, Vc, nanmean(Vc,2));
    Vc = reshape(Vc, nDims, B, C);

    % get variance map for trial aveage
    segFrames = cumsum(floor(segIdx * opts.frameRate)); %max nr of frames per segment
    trialAvg = nanmean(rateDisc_getBhvRealignment(Vc, bhv, segFrames, opts), 3); %align to segments and average over all trials
    covVc = cov(trialAvg(:, ~isnan(trialAvg(1,:)))');  % S x S
    U1 = arrayShrink(U(1:size(allenMask,1), 1:size(allenMask,2), 1:nDims), allenMask, 'merge');
    avgVarMap = nansum((U1 * covVc) .* U1, 2)';  % explained variance map
    
    % get variance map for all data
    Vc = reshape(Vc(1:nDims,:,:), nDims,[]);
    covVc = cov(Vc(:, ~isnan(Vc(1,:)))');  % S x S
    varMap = nansum((U1 * covVc) .* U1, 2)';  % explained variance map
    
    % go through PCs and compute explained variance
    U = U(yRange, xRange, :); %crop to mask range
    U = arrayShrink(U, mask, 'merge');
    
    compVar = zeros(1, size(Vc,1));
    for x = 1 : size(Vc,1)
        cVar = nansum((U(:,1:x) * covVc(1:x,1:x)) .* U(:,1:x), 2);  % explained variance for current nr of components
        compVar(x) = nansum(cVar, 1);
    end
    save([fPath 'compVar.mat'],'compVar','varMap','avgVarMap'); %save explained variance and final variance maps

    sessionVar(:, iRecs) = compVar; %keep explained variance for current session
    sessionVarMap(:, iRecs) = varMap; %keep explained variance map for current session
    sessionAvgVarMap(:, iRecs) = avgVarMap; %keep explained variance map for current session
    fprintf('Recording %d/%d finished.\n', iRecs, length(recs));
    toc
end

%% save output
save([bPath 'sessionVar' trainingRange '.mat'], 'sessionVar', 'sessionVarMap', 'sessionAvgVarMap'); %save explained variance per Vc component for all recordings
