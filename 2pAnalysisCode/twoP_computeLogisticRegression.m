function twoP_computeLogisticRegression(animal,session,num_reps,epochs_only,sub_red,useTrials)
%% Logistic regression - prcesses single session 

% animal = 'CSP30'; session = '200301a';
animal = 'Fez57'; session = '20200716a';
% epochs_only = 1;
% sub_red = 5;
% num_reps = 50;

addpath(genpath(pwd));

S = twoP_settings;

if ~exist('num_reps','var') || isempty(num_reps)
    num_reps = 20;
end

nRep = num_reps;

if ~exist('useTrials','var') || isempty(useTrials)
    useTrials = 250;
end

useTrials = twoP_computeTrialCounts(animal,session,useTrials);

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

if ~exist('epochs_only','var')
    cVc = Vc;
elseif isempty(epochs_only)
    cVc = Vc;
elseif epochs_only > 0
    eVc = twoP_epochTrialMean(Vc,S.allEpoches.idx);
    cVc = eVc;
end

% Define some input parameters

clear lr;

regType = 'lasso'; stepSize = []; decType = 'allChoice'; 
learnType = 'logistic';

%% All neurons, to compute shuffle, multiple iterations will be needed
disp('Running logistic regression choice decoder using entire neural population...')
tic
[lr.cvAcc, lr.bMaps, lr.mdlAll, lr.trialCnt, lr.cvAccShuf] = ...
    rateDisc_logDecoder(cVc, [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,20);
disp(['Running logistic regression using all neurons of this session... done for: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Single neurons
disp('Running logistic regression choice decoder using individual neurons...')

tic
cvAccSingle = zeros(size(cVc,1),size(cVc,2));
parfor i = 1:sum(idxCell(:,1))
    disp(['Running logistic regression for neuron ' num2str(i) '.'])
[cvAccSingle(i,:), ~, ~, ~, ~] = ...
    rateDisc_logDecoder(cVc(i,:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccSingle = cvAccSingle; clear cvAccSingle
disp(['Running logistic regression on individual neurons... done for: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Red (tdT+) neuronal population
disp('Running logistic regression choice decoder using the tdT+ population...')

tic
[lr.cvAccRed, lr.bMapsRed, ~, ~, lr.cvAccRedShuf] = ...
    rateDisc_logDecoder(cVc(idxRed,:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,20);
disp(['Running logistic regression on tdT+ population... done for: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Sequentially adding red (tdT+) neurons (based on their brightness) to the tdT+ population

disp('Running logistic regression on accumulating tdT+ neural population...')
tic
cvAccRedAccum = zeros(length(idxRedSorted),size(cVc,2));
parfor i = 1:length(idxRedSorted)
[cvAccRedAccum(i,:), ~, ~, ~, ~] = ...
    rateDisc_logDecoder(cVc(idxRedSorted(1:i,1),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,0);
end
lr.cvAccRedAccum = cvAccRedAccum;
disp(['Running logistic regression on accumulating tdT+ neural population... done for: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

%% Compute the predictive accuracy of the logistic regression model for
% Number of non-red cells matched to an equal number of red cells

if sub_red > 0
    numRed = sub_red;

    tic
    idxR=find(idxRed);
    idxAll = arrayfun(@(x) randperm(x,numRed),ones(1,nRep)*length(idxR),'UniformOutput',false);
    idxRnew = cell2mat(cellfun(@(x) idxR(x),idxAll,'UniformOutput',false));
    cvAccR = zeros(nRep,size(cVc,2));
    parfor i = 1:nRep
        disp(['Iteration #' num2str(i)])
        [cvAccR(i,:), ~, ~, ~, cvAccRShuf(i,:)] = ...
            rateDisc_logDecoder(cVc(idxRnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,1);
    end
    lr.cvAccR = cvAccR; lr.cvAccRShuf = cvAccRShuf; clear cvAccR
    disp(['Logistic regression completed for tdT+ (subsampled) neurons: ' animal ' ' session  ' in ' num2str(toc) ' seconds.']);

elseif ~exist('sub_red','var') || isempty(sub_red)
    numRed = sum(idxRed);
    idxRnew = repmat(find(idxRed),1,nRep);
end

idxU=find(~idxRed);

tic
idxAll = arrayfun(@(x) randperm(x,numRed),ones(1,nRep)*length(idxU),'UniformOutput',false);
idxUnew = cell2mat(cellfun(@(x) idxU(x),idxAll,'UniformOutput',false));
cvAccU = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    disp(['Iteration #' num2str(i)])
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
    disp(['Iteration #' num2str(i)])
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
    disp(['Iteration #' num2str(i)])
    [cvAccMixedUU(i,:), ~, ~, ~, cvAccMixedUUShuf(i,:)] = ...
        rateDisc_logDecoder(cVc(idxUUnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,learnType,1);
end
lr.cvAccMixedUU = cvAccMixedUU; lr.cvAccMixedUUShuf = cvAccMixedUUShuf; clear cvAccMixedUU
disp(['Logistic regression completed for unlabeled (mixed, shuffled) neurons: ' animal ' ' session ' in ' num2str(toc) ' seconds.']);

lr.useTrials = useTrials;

saveDir= fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir);
if useTrials == 250
    saveFileName = ['LR_' num2str(size(lr.cvAcc,2)) '.mat'];
else 
    saveFileName = ['LR_' num2str(size(lr.cvAcc,2)) '_useTrial' num2str(useTrials) '.mat'];
end
 
save(fullfile(saveDir,saveFileName),'lr');
disp(['Logistic regression data saved as: ' fullfile(saveDir,saveFileName)]);