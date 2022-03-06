function twoP_computeLogisticRegression(animal,session,num_reps,epochs_only,sub_red)
%% 
% animal = 'CSP30'; session = '200301a';
% animal = 'CSP30'; session = '200320';
% epochs_only = 1;
% sub_red = 5;
% num_reps = 50;

addpath(genpath(pwd));

S = twoP_settings;

if ~exist('num_reps','var') || isempty(num_reps)
    num_reps = 20;
end
nRep = num_reps;

imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;

% Load event-aligned imaging data, behavior data, cell-type ID data
Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mat'),'Vc'); Vc = Vc.Vc;
cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
idxCell = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'iscell.npy'));
idxRed = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'redcell.npy'));
idxRed = idxRed(logical(idxCell(:,1)),:);
idxRedSorted = [find(idxRed(:,1)) idxRed(idxRed(:,1)==1,2)]; idxRedSorted = sortrows(idxRedSorted,2);
idxRed = logical(idxRed(:,1));

if exist('epochs_only','var') && (epochs_only > 0)
    eVc = twoP_epochTrialMean(Vc,S.allEpoches.idx);
    cVc = eVc;
else 
    cVc = Vc;
end

%% Logistic regression

% Define some input parameters

clear lr;

useTrials = 250; 
useTrials = twoP_computeTrialCounts(animal,session,useTrials);
regType = 'lasso'; stepSize = []; decType = 'allChoice'; 
learnType = 'logistic';

%% All neurons, to compute shuffle, multiple iterations will be needed
tic
[lr.cvAcc, lr.bMaps, lr.mdlAll, lr.trialCnt, lr.cvAccShuf] = ...
    rateDisc_logDecoder(cVc, [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,20);
disp(['Logistic regression completed for all neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Single neurons 
tic
cvAccSingle = zeros(size(cVc,1),size(cVc,2));
parfor i = 1:length(idxRedSorted)
[cvAccSingle(i,:), ~, ~, ~, ~] = ...
    rateDisc_logDecoder(cVc(i,:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccSingle = cvAccSingle; clear cvAccSingle
disp(['Logistic regression completed for tdT+ neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Red (tdT+) neurons - all
tic
[lr.cvAccRed, lr.bMapsRed, ~, ~, lr.cvAccRedShuf] = ...
    rateDisc_logDecoder(cVc(idxRed,:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,20);
disp(['Logistic regression completed for tdT+ neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Sequentially adding red (tdT+) neurons based on their brightneess
tic
cvAccRedAccum = zeros(length(idxRedSorted),size(cVc,2));
parfor i = 1:length(idxRedSorted)
[cvAccRedAccum(i,:), ~, ~, ~, ~] = ...
    rateDisc_logDecoder(cVc(idxRedSorted(1:i,1),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccRedAccum = cvAccRedAccum;
disp(['Logistic regression completed for tdT+ neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);


%% Compute the predictive accuracy of the logistic regression model for
% non-red cells matched to the number of red cells

if exist('sub_red','var') && (sub_red > 0)
    numRed = sub_red;

    tic
    idxR=find(idxRed);
    idxAll = arrayfun(@(x) randperm(x,numRed),ones(1,nRep)*length(idxR),'UniformOutput',false);
    idxRnew = cell2mat(cellfun(@(x) idxR(x),idxAll,'UniformOutput',false));
    cvAccR = zeros(nRep,size(cVc,2));
    parfor i = 1:nRep
        [cvAccR(i,:), ~, ~, ~, cvAccRShuf(i,:)] = ...
            rateDisc_logDecoder(cVc(idxRnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,1);
    end
    lr.cvAccR = cvAccR; lr.cvAccRShuf = cvAccRShuf; clear cvAccR
    disp(['Logistic regression completed for tdT+ (subsampled) neurons: ' animal ' ' session  ' in ' num2str(toc) ' seconds.']);

else 
    numRed = sum(idxRed);
    idxRnew = repmat(find(idxRed),1,nRep);
end

idxU=find(~idxRed);

tic
idxAll = arrayfun(@(x) randperm(x,numRed),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUnew = cell2mat(cellfun(@(x) idxU(x),idxAll,'UniformOutput',false));
cvAccU = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccU(i,:), ~, ~, ~, cvAccUShuf(i,:)] = ...
        rateDisc_logDecoder(cVc(idxUnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,1);
end
lr.cvAccU = cvAccU; lr.cvAccUShuf = cvAccUShuf; clear cvAccU
disp(['Logistic regression completed for unlabeled (population-matched) neurons: ' animal ' ' session  ' in ' num2str(toc) ' seconds.']);

tic
idxAll = arrayfun(@(x) randperm(x,numRed),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUnew = cell2mat(cellfun(@(x) idxU(x),idxAll,'UniformOutput',false));
cvAccMixedUR = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccMixedUR(i,:), ~, ~, ~, cvAccMixedURShuf(i,:)] = ...
        rateDisc_logDecoder(cVc([idxUnew(:,i);idxRnew(:,i)],:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,1);
end
lr.cvAccMixedUR = cvAccMixedUR; lr.cvAccMixedURShuf = cvAccMixedURShuf; clear cvAccMixedUR
disp(['Logistic regression completed for mixed (shuffled) neurons: ' animal ' ' session  ' in ' num2str(toc) ' seconds.']);

tic
idxUU = arrayfun(@(x) randperm(x,2*numRed),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUUnew = cell2mat(cellfun(@(x) idxU(x),idxUU,'UniformOutput',false));
cvAccMixedUU = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccMixedUU(i,:), ~, ~, ~, cvAccMixedUUShuf(i,:)] = ...
        rateDisc_logDecoder(cVc(idxUUnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,1);
end
lr.cvAccMixedUU = cvAccMixedUU; lr.cvAccMixedUUShuf = cvAccMixedUUShuf; clear cvAccMixedUU
disp(['Logistic regression completed for unlabeled (mixed, shuffled) neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

lr.useTrials = useTrials;

saveDir= fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir);
saveFileName = ['LR_' num2str(size(lr.cvAcc,2)) '.mat'];
save(fullfile(saveDir,saveFileName),'lr');
disp(['Logistic regression data saved as: ' fullfile(saveDir,saveFileName)]);