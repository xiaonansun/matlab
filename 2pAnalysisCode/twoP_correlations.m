animal = 'Plex50'; session = '200330a';

S = twoP_settings;
imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;

sPerFrame = S.msPerFrame/1000;

% Load event-aligned imaging data, behavior data, cell-type ID data
Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mat'),'Vc'); Vc = Vc.Vc;
cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
idxCell = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'iscell.npy'));
idxRed = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'redcell.npy'));
idxRed = logical(idxRed(logical(idxCell(:,1))));

% Define the epoches
idxEpochInit = [S.segFrames(1) S.segFrames(2)];
idxEpochStim = [S.segFrames(2)+1 S.segFrames(3)];
idxEpochStimEarly = [S.segFrames(2)+1 floor(mean(idxEpochStim))];
idxEpochStimLate = [ceil(mean(idxEpochStim)) S.segFrames(3)];
idxEpochDelay = [S.segFrames(3)+1 S.segFrames(4)];
idxEpochResponse = [S.segFrames(4)+1 S.segFrames(5)];
idxEpochResponseEarly = [S.segFrames(4)+1 floor(mean(idxEpochResponse))];
idxEpochResponseLate = [ceil(mean(idxEpochResponse)) S.segFrames(5)];
idxEpochAll = [idxEpochInit; idxEpochStimEarly; idxEpochStimLate; idxEpochDelay; idxEpochResponseEarly; idxEpochResponseLate];

% Extract the indices of specific type of trials (i.e all right-sided response trials)
bhv = twoP_bhvSubSelection(cBhv);
stimEL = squeeze(mean(Vc(:,idxEpochStim(1):idxEpochStim(2),bhv.stim.EasyLeft),2,'omitnan'));
stimER = squeeze(mean(Vc(:,idxEpochStim(1):idxEpochStim(2),bhv.stim.EasyRight),2,'omitnan'));
stimAL = squeeze(mean(Vc(:,idxEpochStim(1):idxEpochStim(2),bhv.stim.AllLeft),2,'omitnan'));
stimAR = squeeze(mean(Vc(:,idxEpochStim(1):idxEpochStim(2),bhv.stim.AllRight),2,'omitnan'));
delayL = squeeze(mean(Vc(:,idxEpochDelay(1):idxEpochDelay(2),bhv.response.Left),2,'omitnan'));
delayR = squeeze(mean(Vc(:,idxEpochDelay(1):idxEpochDelay(2),bhv.response.Right),2,'omitnan'));
responseL = squeeze(mean(Vc(:,idxEpochResponse(1):idxEpochResponse(2),bhv.response.Left),2,'omitnan'));
responseR = squeeze(mean(Vc(:,idxEpochResponse(1):idxEpochResponse(2),bhv.response.Right),2,'omitnan'));

%% 

eVc = twoP_epochTrialMean(Vc,idxEpochAll);

%% Logistic regression
cVc = eVc;
clear lr;
useTrials = 250; regType = 'lasso'; stepSize = []; decType = 'allChoice';

% All neurons
[lr.cvAcc, lr.bMaps, lr.mdlAll, lr.trialCnt, lr.cvAccShuf, lr.bMapsShuf, lr.mdlAllShuf, lr.leftIdxShuf, lr.leftIdx] = rateDisc_logDecoder(cVc, [], cBhv, useTrials, 0, regType, stepSize, decType,[]);

% Red neurons
[lr.cvAccRed, lr.bMapsRed, ~, ~, ~, ~, ~, ~, ~] = rateDisc_logDecoder(cVc(idxRed,:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,[]);

% Compute the predictive accuracy of the logistic regression model for
% non-red cells matched to the number of red cells
clear idxNRnew cvAccNR
nRep = 50;
idxNR=find(~idxRed);
% idxNRnew = zeros(sum(idxRed),nRep);
idxAll = arrayfun(@(x) randperm(x,sum(idxRed)),ones(1,nRep)*length(idxNR),'UniformOutput',false);
idxNRnew = cell2mat(cellfun(@(x) idxNR(x),idxAll,'UniformOutput',false));
cvAccNR = zeros(nRep,size(cVc,2));
parfor i = 1:nRep
    [cvAccNR(i,:), ~, ~, ~, ~, ~, ~, ~, ~] = rateDisc_logDecoder(cVc(idxNRnew(:,i),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType,[]);
end
lr.cvAccNR = cvAccNR; clear cvAccNR


%%
eVc = twoP_epochTrialMean(Vc,idxEpochAll);
[lr.cvAcc, lr.bMaps, lr.mdlAll, lr.trialCnt, lr.cvAccShuf, lr.bMapsShuf, lr.mdlAllShuf, lr.leftIdxShuf, lr.leftIdx] = rateDisc_logDecoder(eVc, [], cBhv, useTrials, 0, regType, stepSize, decType,[]);

%% Logistic regression - predictive accuracy of individual neurons
singleAcc = zeros(size(eVc,1),size(eVc,2));
parfor i = 1:size(singleAcc,1)
%     defVc = eVc(i,:,:); 
defVc = eVc;
    defVc(i,:,:)=[];
    [singleAcc(i,:),~,~,~,~,~,~,~] = rateDisc_logDecoder(defVc, [], cBhv, useTrials, 0, regType, stepSize, decType,[]);
    disp(num2str(i));
end
diffAcc = lr.cvAcc - singleAcc;
%%
errorShuf = std(lr.shufAUC,0,3,'omitnan');
plot(lr.allAUC'); hold on;
plot(0.5+2*errorShuf','k'); 
plot(0.5-2*errorShuf','k'); 

%% AUC analysis
% eVc = zeros(size(Vc,1),size(idxEpochAll,1),size(Vc,3));
% for i = 1:size(eVc,2)
%     for j = 1:size(Vc,3)
%         eVc(:,i,j)= mean(Vc(:,idxEpochAll(i,1):idxEpochAll(i,2),j),2,'omitnan');
%     end
% end

eVc = twoP_epochTrialMean(Vc,idxEpochAll);
[auc.all, auc.shuffle, auc.sAUC] = twoP_aucAnalysisNew(eVc, cBhv, 1, 1, animal, session);

%% 
numBins = 30;
limX = [0 1];
sigAUC = nan(size(auc.all,1),size(auc.all,2));
binCounts = zeros(size(auc.all,2),numBins);
uShuf = mean(auc.shuffle,3)+2*std(auc.shuffle,0,3);
lShuf = mean(auc.shuffle,3)-2*std(auc.shuffle,0,3);
uAUC = (auc.all >= uShuf).*auc.all;
lAUC = (auc.all <= lShuf).*auc.all;
sigAUC = uAUC + lAUC; sigAUC(sigAUC==0) = nan;

figure('Position',[500 500 1200 250])
tiledlayout(1,size(auc.all,2))
for i = 1:size(sigAUC,2)
    binCounts(i,:) = histcounts(sigAUC(:,i), 30,'BinLimits',limX);
end
limY = [min(binCounts(:)) ceil(max(binCounts(:))/10)*10];

for i = 1:size(sigAUC,2)
    nexttile
    histogram(sigAUC(~isnan(sigAUC(:,i)),i), numBins,'BinLimits',limX);
    ylim([limY(1) limY(end)]);
end

%% Plot AUC analysis

mEpochAUC = zeros(size(auc.all,1),size(idxEpochAll,1));
mEpochAUCshuf = zeros(size(auc.all,1),size(idxEpochAll,1),size(auc.shuffle,3));

for i = 1:size(idxEpochAll,1)
    mEpochAUC(:,i)= mean(auc.all(:,idxEpochAll(i,1):idxEpochAll(i,2)),2);
    for j = 1:size(auc.shuffle,3)
        mEpochAUCshuf(:,i,j) = mean(auc.shuffle(:,idxEpochAll(i,1):idxEpochAll(i,2),j),2);
    end
end

mEpochAUCshufMean = mean(mEpochAUCshuf,3);
mEpochAUCshufStd = std(mEpochAUCshuf,0,3);



%%
[svm.cvAcc, svm.bMaps, svm.betaNeuron, svm.mdlAll, svm.trialCnt, svm.allAUC, svm.shufAUC] = rateDisc_logDecoder(Vc, [], cBhv, useTrials, 0, regType, stepSize, decType,'svm',1);

lr.sorted.betaDelay = mean(lr.bMaps(:,idxEpochDelay(1):idxEpochDelay(2)),2,'omitnan');
[lr.sorted.betaDelay,lr.sorted.ibetaDelay]=sort(lr.sorted.betaDelay);
svm.sorted.betaDelay = mean(svm.bMaps(:,idxEpochDelay(1):idxEpochDelay(2)),2,'omitnan');
[svm.sorted.betaDelay,svm.sorted.ibetaDelay]=sort(svm.sorted.betaDelay);

auc.sorted.delay = mean(lr.allAUC(:,idxEpochDelay(1):idxEpochDelay(2)),2,'omitnan');
[auc.sorted.delay,auc.sorted.idelay]=sort(auc.sorted.delay);
imagesc(squeeze(Vc(auc.sorted.idelay(end),:,:))');

%%
pos = lr.allAUC > (mean(lr.shufAUC,3) + std(lr.shufAUC,0,3)*2);
neg = lr.allAUC < (mean(lr.shufAUC,3) - std(lr.shufAUC,0,3)*2);

posRed = sum(pos(idxRed,:)); negRed = sum(neg(idxRed,:));
plot(posRed,'k'); hold on; plot(negRed,'b');

figure(2)
posNonred = sum(pos(~idxRed,:)); negNonred = sum(neg(~idxRed,:));
plot(posNonred,'k'); hold on; plot(negNonred,'b');
%%
plot(mean(squeeze(Vc(auc.sorted.idelay(round(end/2)),:,bhv.response.Right))',1,'omitnan'))
hold on;
plot(mean(squeeze(Vc(auc.sorted.idelay(round(end/2)),:,bhv.response.Left))',1,'omitnan'))
%%
stimInput = responseL;
% stimEL = [stimEL((idxRed),:);stimEL(~idxRed,:)]; % sort by red-red,
% red-nonred and nonred-nonred correlations

% [R,P] = corrcoef(stimEL');
for i = 1:size(stimInput,1)
    stimInput(i,isoutlier(stimInput(i,:)))=nan; % Remove outliers which would otherwise falsely elevate correlation coefficients
end

[R,P] = corrcoef(stimInput','Rows','pairwise');
idxCol = repmat((1:length(R))',length(R),1);
idxRow = repmat(1:length(R),length(R),1); idxRow = idxRow(:);
A = tril(R,-1); A(A==0)=nan;
[B,iB] = sort(A(:),'descend'); iB(isnan(B))=[]; B(isnan(B))=[];


RR = A(idxRed,:);
RN = A(~idxRed,idxRed);
NN = A(~idxRed,~idxRed);cdfplot(RR(:))
hold on;
cdfplot(RN(:))
cdfplot(NN(:))

% fig_configAxis

[h,p] = kstest2(NN(:),RR(:))
[h,p] = kstest2(NN(:),RN(:))
[h,p] = kstest2(RR(:),RN(:))
legend({'Pos-Pos';'Pos-Neg';'Neg-Neg'})


%%
imagesc(R)
imagesc(P)