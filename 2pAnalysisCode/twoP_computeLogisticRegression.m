function twoP_computeLogisticRegression(animal,session,num_reps)
%% 
% animal = 'Plex50'; session = '200330a';


if ~exist('num_reps','var') || isempty(num_reps)
    num_reps = 20;
end
nRep = num_reps;


S = twoP_settings;
imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;


% Load event-aligned imaging data, behavior data, cell-type ID data
Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mat'),'Vc'); Vc = Vc.Vc;
cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
idxCell = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'iscell.npy'));
idxRed = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'redcell.npy'));
idxRed = logical(idxRed(logical(idxCell(:,1))));

%% Logistic regression
% cVc = eVc;
cVc = Vc;

clear lr;
useTrials = 250; regType = 'lasso'; stepSize = []; decType = 'allChoice'; 
learnType = 'logistic';

% All neurons, to compute shuffle, multiple iterations will be needed
tic
[lr.cvAcc, lr.bMaps, lr.mdlAll, lr.trialCnt, lr.cvAccShuf, lr.bMapsShuf] = ...
    rateDisc_logDecoder(cVc, [], cBhv, useTrials, 0, regType, stepSize, decType,learnType);
disp(['Logistic regression completed for all neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

% Red neurons
tic
[lr.cvAccRed, lr.bMapsRed, ~, ~, ~, ~] = ...
    rateDisc_logDecoder(cVc(idxRed,:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
disp(['Logistic regression completed for tdT+ neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

% Compute the predictive accuracy of the logistic regression model for
% non-red cells matched to the number of red cells

tic
idxU=find(~idxRed);
idxAll = arrayfun(@(x) randperm(x,sum(idxRed)),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUnew = cell2mat(cellfun(@(x) idxU(x),idxAll,'UniformOutput',false));
cvAccU = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccU(i,:), ~, ~, ~, ~, ~] = ...
        rateDisc_logDecoder(cVc(idxUnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccU = cvAccU; clear cvAccU
disp(['Logistic regression completed for unlabeled (population-matched) neurons: ' animal ' ' session  ' in ' num2str(toc) ' seconds.']);

tic
idxAll = arrayfun(@(x) randperm(x,sum(idxRed)),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUnew = cell2mat(cellfun(@(x) idxU(x),idxAll,'UniformOutput',false));
cvAccMixedUR = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccMixedUR(i,:), ~, ~, ~, ~, ~] = ...
        rateDisc_logDecoder(cVc([idxUnew(:,i);find(idxRed)],:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccMixedUR = cvAccMixedUR; clear cvAccMixedUR
disp(['Logistic regression completed for mixed (shuffled) neurons: ' animal ' ' session]);

tic
idxUU = arrayfun(@(x) randperm(x,2*sum(idxRed)),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUUnew = cell2mat(cellfun(@(x) idxU(x),idxUU,'UniformOutput',false));
cvAccMixedUU = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccMixedUU(i,:), ~, ~, ~, ~, ~] = ...
        rateDisc_logDecoder(cVc(idxUUnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccMixedUU = cvAccMixedUU; clear cvAccMixedUU
disp(['Logistic regression completed for unlabeled (mixed, shuffled) neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

saveDir= fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir);
saveFileName = ['LR_' num2str(size(lr.cvAcc,2)) '.mat'];
save(fullfile(saveDir,saveFileName),'lr');
disp(['Logistic regression data saved as: ' fullfile(saveDir,saveFileName)]);