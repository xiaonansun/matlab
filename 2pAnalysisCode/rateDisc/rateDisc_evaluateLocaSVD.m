function rateDisc_evaluateLocaSVD(animal,nDims)
%% some variables
% nDims = [100, 500, 1000]; %number of whole-frame components to test
if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
dPath = [cPath animal filesep 'SpatialDisc' filesep]; % path for raw data

stepSize = 10000;
frameOffset = 50; %skip the first couple of frames to avoid filter artefact

%% block info
load([bPath 'trialInfo.mat'], 'trialCnt', 'recs');

load([bPath 'locAC.mat'],'A','C','nanIdx');
nanIdx = reshape(nanIdx, [], sum(trialCnt));

load([bPath 'mask.mat'],'allenMask','xRange','yRange');
mask = isnan(nanmean(A,3)) | allenMask(yRange,xRange);
A = arrayShrink(A,mask,'merge');

%sort A by highest eigenvalues
load([bPath 'locanmf_decomp_loc50.mat'],'lambdas');
[~,b] = sort(lambdas,'ascend');
A = A(:,b);
C = C(b,:);

%% go through recordings and compute variance and prediction error
allMse = NaN(length(nDims),length(recs));
for iRecs = 2 : length(recs)
    tic
    fPath = [dPath filesep recs(iRecs).name filesep]; disp(fPath); %session data path
    load([fPath 'Vc.mat'],'Vc','U');
    load([fPath 'opts2.mat'],'opts');
    
    U = alignAllenTransIm(single(U),opts.transParams); %align to allen
    U = U(yRange, xRange, :); %crop to mask range
    
    cIdx = sum(trialCnt(1:iRecs-1))+1 : sum(trialCnt(1:iRecs)); %trials for current recording
    
    prevFrames = ~nanIdx(:, 1 : cIdx(1)-1); %all frames until current trials
    currFrames = ~nanIdx(:, cIdx); %all frames until current trials
    wV = NaN([size(C,1) size(currFrames)], 'single');
    wV(:,currFrames) = C(:,sum(prevFrames(:)) + 1 : sum(prevFrames(:)) + sum(currFrames(:))); %locaNMF for current session
    
    Vc = reshape(Vc(:,~isnan(Vc(1,:))),size(Vc,1),[]); %exclude NaN frames
    wV = reshape(wV(:,~isnan(wV(1,:))),size(wV,1),[]); %exclude NaN frames
    
    wV = bsxfun(@minus, wV, mean(wV,2)); %make sure wV is zero-mean
    Vc = bsxfun(@minus, Vc, mean(Vc,2)); %make sure Vc is zero-mean
    
    for iDims = 1 : length(nDims)
        if nDims(iDims) <= size(A,2)
            for iFrames = frameOffset : stepSize : size(wV,2)
                
                if size(wV,2) > iFrames + stepSize
                    cData = A(:,1:nDims(iDims)) * wV(1:nDims(iDims),iFrames+1 : iFrames + stepSize);
                else
                    cData = A(:,1:nDims(iDims)) * wV(1:nDims(iDims),iFrames+1 : end);
                end
                
                % first run for current recording
                if ndims(U) == 3
                    U = arrayShrink(U, mask, 'merge');
                    varMap = zeros(size(U,1),1);
                    errMap = zeros(size(U,1),1);
                end
                if size(Vc,2) > iFrames + stepSize
                    rData = U * Vc(:, iFrames+1 : iFrames + stepSize); %reconstruct 'raw' data
                else
                    rData = U * Vc(:, iFrames+1 : end); %reconstruct 'raw' data
                end
                
                varMap = varMap + sum(rData.^2,2); %sum variance in raw data
                errMap = errMap + sum((rData-cData).^2,2); %sum squared prediction error
                clear rData cData
            end
            
            varMap = varMap ./ (size(wV,2) - frameOffset); %variance for each pixel in raw data
            errMap = errMap ./ (size(wV,2) - frameOffset); %squared prediction error in each pixel
            
            cVar = sum(varMap); %total variance
            errVar = sum(errMap); %total squared error
            save([fPath 'locaBlockVar_'  num2str(nDims(iDims)) '.mat'],'cVar','varMap','errVar','errMap');
            
            allMse(iDims,iRecs) = errVar / cVar; %percent error
            fprintf('Used dims: %d; loss: %f percent.\n', nDims(iDims), allMse(iDims,iRecs)*100);
        end
    end
    fprintf('Recording %d/%d finished.\n', iRecs, length(recs));
    toc
end

%% save output
for iDims = 1 : length(nDims)
    mse = allMse(iDims,:);
    save([bPath 'locaMse_'  num2str(nDims(iDims)) '.mat'], 'mse');
end
