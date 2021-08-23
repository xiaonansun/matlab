% function twoP_scripts(animal,session)
% mfilename('fullpath')
%% Load trial-aligned fluorescence data
clear animal session data npy SessionData events; close all;

% Specify session to load
animal = 'CSP27';
session = '20200321a';
[npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData(animal,session, true, 10, false);

% baseFileName = [animal '_' session];
% ----- Makes adjustments to the twoP data struct ----- %%
data = twoP_adjustData(data,SessionData);

% ----- Define event-aligned matrices ----- %%
[lick,data.lickWinIdx,data.lickWinMs,data.dataLick, data.dataLickTrialNumbers]=twoP_alignToLick(data, SessionData);

% ----- Behavior analysis----- %%
% events =SessionData.stimEvents;
% % organize stimulus events of each trial (row) into columns: column 1 is
% % left and column 2 is right
% events = cellfun(@(x) x{:},[cellfun(@(x) x(1),events,'UniformOutput',false)' cellfun(@(x) x(2),events,'UniformOutput',false)'],'UniformOutput',false); % In this line, we count the number of auditory events emerging from the left or right speakers and report them in two columns (column 1 for left and 2 for right)
% frEvents= cellfun(@numel,events);

E = behavior_getStimEvents(SessionData); events = E.events; frEvents = E.frEvents; ratioEvents = E.ratioEvents;
% frEvents(:,4)=(SessionData.CorrectSide)';
% frEvents(:,5)=(SessionData.DistStim)'; % distractor

%% Logistic regression decoder 

% [cvAcc, bMaps, trialCnt] = rateDisc_logDecoder(data, U, bhv, useTrials, targMod, regType, stepSize, decType)
% useTrials: 400 trials is default, can be adjusted, fewer trials will reduce the accuracy of the decoder
% stepSize: will downsample data
% decType: the type of decoder, default is allchoice
% Vc: neural data neurons x frame x trial
% cBhv: SessionData (Behai vior data aka bpod output)
% segFrames: defined at the top of this cell
% opts: struct defined at the top of this cell

% dbstop if error

lr = twoP_logisticRegression(animal, session, data, SessionData, true, true); %% the last variable, if true, loads existing logstic regression data

cBhv = selectBehaviorTrials(SessionData, data.trialNumbers); %% very important in matching trial indices

%% 

close all;

Vc = lr.Vc; Vc_nr = lr.Vc_nr; Vc_r = lr.Vc_r;
peth.mVc = twoP_sortROIbyActivity(mean(Vc,3,'omitnan'),[],'sgolay',1);
peth.mVc_nr = twoP_sortROIbyActivity(mean(Vc_nr,3,'omitnan'),[],'sgolay',1);
peth.mVc_r = twoP_sortROIbyActivity(mean(Vc_r,3,'omitnan'),[],'sgolay',1);

fig = figure(1);
fig.Position = [100,100,400,800];

tiledlayout(4,1,...
    'TileSpacing','none');
nexttile([3 1]); 
imagesc([peth.mVc_r; peth.mVc_nr]); hold on;
yline(size(peth.mVc_r,1),'-w');
xline(lr.segFrames+1,'-w');
ax = fig.CurrentAxes; ax.XTick = [];

nexttile
plot(mean(peth.mVc,1),'-k'); hold on;
plot(mean(peth.mVc_nr,1),'-g');
plot(mean(peth.mVc_r,1),'-r');

xlim([1 size(Vc,2)]);
xline(lr.segFrames+1,'-k');
ax = fig.CurrentAxes;
ax = fig_configAxis(ax);
% figure(2);
% plot(mean(twoP_sortROIbyActivity(mean(Vc,3))));
% imagesc(mean(Vc,3))

%% Correlations (noise and signal)

easyL = squeeze(sum(data.neural(:,data.trialStimFrame:data.trialStimFrame+round(1000/data.msPerFrame),ratioEvents==0),2));
easyR = squeeze(sum(data.neural(:,data.trialStimFrame:data.trialStimFrame+round(1000/data.msPerFrame),ratioEvents==1),2));
easyLr = easyL(data.idx_redcell,:);
easyLnr = easyL(data.idx_notredcell,:);
easyLsorted = [easyLr;easyLnr];
% easyL = easyLr; clear easyLr;

R_L=corrcoef(easyL');
R_Lr = corrcoef(easyLr');
R_Lnr = corrcoef(easyLnr');
R_Lsorted=corrcoef(easyLsorted');

% Red cells
R_Lr_upper = triu(R_Lr,1); R_Lr_upper(R_Lr_upper ==0) = NaN;
R_Lr_upper_vec = reshape(R_Lr_upper,[size(R_Lr,1)*size(R_Lr_upper,2),1]);
R_Lr_upper_vec(isnan(R_Lr_upper))=[];

% Non-red cells
R_Lnr_upper = triu(R_Lnr,1); R_Lnr_upper(R_Lnr_upper ==0) = NaN;
R_Lnr_upper_vec = reshape(R_Lnr_upper,[size(R_Lnr,1)*size(R_Lnr_upper,2),1]);
R_Lnr_upper_vec(isnan(R_Lnr_upper))=[];

% Calculate distance between cells
M_dist = zeros(size(R_L,1),size(R_L,2));
for i = 1:size(R_L,1)
    for j = 1:size(R_L,2)
        M_dist(i,j)=pdist2(npy.stat{i}.med,npy.stat{j}.med);
    end
end

histogram(M_dist(:))

figure(1)
imagesc(R_L);

figure(2)
imagesc(R_Lsorted);

figure(3)
imagesc(R_Lr);

figure(4)
imagesc(R_Lr_upper)

figure(5)
imagesc(R_Lnr_upper)

figure(6)
histogram(R_Lr_upper_vec,'NumBins',30,'Normalization','pdf'); hold on;
histogram(R_Lnr_upper_vec,'Normalization','pdf');

%% Compute basic firing rate statistics
sPF = data.msPerFrame/1000;
tVec = sPF:sPF:sPF*size(npy.spks,2);
tVecRounded = floor(tVec);
for i = 0:max(tVecRounded)-1
    spksPerSec(:,i+1) = sum(npy.spks(:,tVecRounded==i),2);
end

% Compute the number of one-second bins of individual neurons for which
% there is no neural activity (spks = 0)
noActivity = sum(spksPerSec==0,2);
histogram(noActivity(npy.iscell(:,1)==1 & npy.redcell(:,1)==0),'Normalization','pdf');
hold on; 
histogram(noActivity(npy.iscell(:,1)==1 & npy.redcell(:,1)==1),'Normalization','pdf');

spksPerSecNRV = reshape(spksPerSec(npy.iscell(:,1)==1 & npy.redcell(:,1)==0,:), 1,sum(npy.iscell(:,1)==1 & npy.redcell(:,1)==0)*size(spksPerSec,2));
spksPerSecRV = reshape(spksPerSec(npy.iscell(:,1)==1 & npy.redcell(:,1)==1,:), 1,sum(npy.iscell(:,1)==1 & npy.redcell(:,1)==1)*size(spksPerSec,2));
histogram(spksPerSecNRV(spksPerSecNRV>0),500,'Normalization','cdf');
hold on;
h= histogram(spksPerSecRV(spksPerSecRV>0),500,'Normalization','cdf');

% durSess = size(npy.spks,2)*data.msPerFrame/1000;
% iFrame = (1:1:durSess)/fps;
%% Reorganize data for PCA and dPCA
% trialID.imaging is the master trial index of 2p imaging data. This
% structure field should contain all indices contained by the other trialID
% fields defined here

% --- Use only one of these --- %
% data.neuralNorm=twoP_normMean(data.neural);
% data.neuralNorm=twoP_normMean(data.neural(data.idx_redcell,:,:));
data.neuralNorm=twoP_normMean(data.neural(data.idx_notredcell(1:length(data.idx_redcell)),:,:));
% --- Use only one of these --- %

data.pcaTime=(data.neuralTimesMs)/1000;
numOfBins = 2;
if numOfBins > 0
    edges = 0:1/numOfBins:1; % divide the auditory stimuli presented throughout the session into bins which in this case is numOfBins
    edges = [0 0.01 edges(2:end-1) 0.99 1];
else 
    edges = [0 0.5 1];
end

stimValues=ratioEvents; % fraction of left sided events

% REMOVE AUTO-REWARD TRIALS
trialID.notAutoReward=find(~SessionData.AutoReward); % Find the IDs of non-auto-reward trials
[trialID.imaging,trialID.trialNumbersIdx,trialID.AutoRewardIdx]= intersect(data.trialNumbers,trialID.notAutoReward); % Find the non-auto-reward imaging trial IDs

% ADDITIONALLY, TRAINING MODE TRIALS ARE ALSO REMOVED: SessionData.SingleSpout
trialID.leftLick = find(SessionData.ResponseSide==1); trialID.rightLick = find(SessionData.ResponseSide==2); % Find the behavior trial IDs of left or right decisions
trialID.imagingLeftLick = intersect(trialID.imaging,trialID.leftLick); % find the trial IDs of non-auto-reward left decisions
trialID.imagingRightLick = intersect(trialID.imaging,trialID.rightLick); % find the trial IDs of non-auto-reward right decisions

[binCounts,edges,stimulusCondition]=histcounts(stimValues(trialID.imaging),edges); % divides the
trialID.stimulusCondition = cell(max(unique(stimulusCondition)),1); % creates a field of the cell class to hold trial ID information for each stimulus category
data.pcaBinCounts = binCounts;

for i = 1:max(unique(stimulusCondition))
    trialID.stimulusCondition{i}=find(stimulusCondition==i);
    trialID.stimulusConditionLeftDec{i}=intersect(trialID.stimulusCondition{i}, trialID.imagingLeftLick);
    trialID.stimulusConditionRightDec{i}=intersect(trialID.stimulusCondition{i}, trialID.imagingRightLick);
end
stimConditions = {trialID.stimulusConditionLeftDec trialID.stimulusConditionRightDec};
% [maxStimTrials,maxStimTrialsIdx]= max(cellfun(@numel ,trialID.stimulusCondition));
[maxRight,maxRightIdx]=max(cellfun(@numel, trialID.stimulusConditionRightDec));
[maxLeft, maxLeftIdx]=max(cellfun(@numel, trialID.stimulusConditionLeftDec));
maxStimTrials = max([maxLeft maxRight]);

% Create an empty multidimensional array for demixed PCA (Kobak et al., 2016): 
% dimensions: (1) neuron ID, (2) time, (3) trials grouped by stimulus
% condition, (4) stimulus condition, (5) decision: left or right
data.pcaArray=nan(size(data.neuralNorm,1),size(data.neuralNorm,2),maxStimTrials,length(trialID.stimulusCondition),2);

% m=0; probably can delete this line
for i = 1:size(data.pcaArray,5) % i: response (left or right)
    for j = 1:size(data.pcaArray,4) % j: stimulus category (1 through 11 with 1=left? and 11=being right? and everything else in between)
        for k = 1:length(stimConditions{i}{j})
            data.pcaArray(:,:,k,j,i) = data.neuralNorm(:,:,ismember(data.trialNumbers,stimConditions{i}{j}(k))); % MODIFY THIS LINE TO INTERCHANGE BETWEEN RED VS NON-RED CELLS
%             m=m+1; probably can delete this line
        end
    end
end

% Permute the dimensions of the PCA array: (1) neuron ID, (2) stimulus conditions, (3) Decision (L/R) grouped by stimulus
% condition, (4) time, (5) trial count
data.pcaArrayPerm=permute(data.pcaArray, [1 4 5 2 3]);
firingRates= data.pcaArrayPerm;
firingRatesAverage = nanmean(firingRates,5);

trialNum=zeros(size(data.pcaArrayPerm,1),size(data.pcaArrayPerm,2),size(data.pcaArrayPerm,3));
for i = 1:size(data.pcaArrayPerm,3)
    for j = 1:size(data.pcaArrayPerm,1)
        for k = 1:size(data.pcaArrayPerm,2)
            trialNum(j,k,i)=sum(sum(~isnan(firingRates(j,k,i,1,:)),5));
        end
    end
end
data.pca.firingRates = firingRates;
data.pca.firingRatesAverage = firingRatesAverage;
data.pca.trialNum=trialNum;

%Run demixed PCA
dpcaTwoPData(data);

%% Re-organize data and plot tuning curves
data.neuralNorm=twoP_normMean(data.neural);
time=data.neuralTimesMs;
stimRange = [500 1000]; % stimulus start and end time (in milliseconds)
idxtime=find(time > stimRange(1) & time <= stimRange(2));
numOfBins = 5; % Number of bins (to be specified by the user) excluding the two extremes: 0 and 1. In total, the number of bins will be numOfBins + 2
edges = 0:1/numOfBins:1; % divide the auditory stimuli presented throughout the session into bins 
edges = [0 0.01 edges(2:end-1) 0.99 1];

stimValues=ratioEvents; % fraction of left sided events
[binCounts,edges,stimulusCondition]=histcounts(stimValues(trialID.imaging),edges); % divides the vector (each entry is a trial) containing the fraction of left-sided events into bins
sR=nan(size(data.neural,1),length(idxtime),max(binCounts),length(binCounts)); % sR = stimulus response. Creates a 4-D vector. 1: cell ID, 2: time, 3: trial, 4: stimulus condition (from all left- to all right-sided events)

for i = 1:length(binCounts)
% sR(:,:,1:binCounts(i),i) = data.neural(:,idxtime,trialID.imaging(stimulusCondition==i));
sR(:,:,1:binCounts(i),i) = data.neural(:,idxtime,stimulusCondition==i);
end

tM = squeeze(nanmean(sum(sR,2),3)); % tM = tuning matrix, where rows are cells and columns are the stimulus conditions from 100% left-sided events to 100% right-sided events.
[sMax,sMaxIdx]=max(tM,[],2);
sMaxIdx(:,2)=1:length(sMaxIdx); sMaxIdx(:,3:4)=sortrows(sMaxIdx(:,1:2),1); %sort the matrix by the stimulus condition with the maximum response
hTuningCurve = figure('Name',['Tuning Curve: ' animal ' ' session]);
imagesc(tM(sMaxIdx(:,4),:)); colormap(jet);
% imagesc(twoP_normMean(tM(sMaxIdx(:,4),:))); colormap(jet);

% No need to remove training mode trials (SessionData.SingleSpout)
% No need to remove auto-reward trials, since we are plotting sensory
% tuning curves
% trialID.notAutoReward=find(~SessionData.AutoReward); % Find the IDs of non-auto-reward trials
% [trialID.imaging,trialID.trialNumbersIdx,trialID.AutoRewardIdx]= intersect(data.trialNumbers,trialID.notAutoReward); % Find the non-auto-reward imaging trial IDs

% trialID.leftLick = find(SessionData.ResponseSide==1); trialID.rightLick = find(SessionData.ResponseSide==2); % Find the behavior trial IDs of left or right decisions
% trialID.imagingLeftLick = intersect(trialID.imaging,trialID.leftLick); % find the trial IDs of non-auto-reward left decisions
% trialID.imagingRightLick = intersect(trialID.imaging,trialID.rightLick); % find the trial IDs of non-auto-reward right decisions

%% Plot stimulus-aligned mean (across trials) traces (PSTH)
sessionDir = suite2pDir(1:strfind(suite2pDir,'suite2p')-1);
sessAvgDir = [sessionDir 'sessAvg' filesep];

if ~exist(sessAvgDir, 'dir')
    mkdir(sessAvgDir)
end

stimWinIdx= data.neuralTimes;
trialID.allImaging=lick(5,:);
trialID.all=intersect(trialID.allImaging,find(~SessionData.AutoReward));
trialID.leftEasy = intersect(trialID.all,find(ratioEvents==0));
[trialID.leftEasyIdx, trialID.leftEasySubIdx,trialID.leftEasyImgIdx]= intersect(trialID.leftEasy,trialID.allImaging);
trialID.rightEasy = intersect(trialID.all,find(ratioEvents==1));
[trialID.rightEasyIdx,trialID.rightEasySubIdx,trialID.rightEasyImgIdx] = intersect(trialID.rightEasy,trialID.allImaging);
correctSide = SessionData.CorrectSide(trialID.all); % correctSide is the side with more stimulus events
sideOneIdx = find(correctSide==1); sideTwoIdx = find(correctSide==2);
% [A,iC,iD]=intersect(data.dataLickTrialNumbers,data.trialNumbers);
stim.unsorted= data.neural;

xTickStim=[1 find(stimWinIdx == 0) find(stimWinIdx == round(1000/data.msPerFrame)) length(stimWinIdx)]; % define the position of x-tick marks
xTickStimLabel = {round(stimWinIdx(1)*data.msPerFrame/1000,2,'significant'), stimWinIdx(stimWinIdx==0), round(stimWinIdx(stimWinIdx == round(1000/data.msPerFrame))*data.msPerFrame/1000,2,'significant'), round(stimWinIdx(end)*data.msPerFrame/1000,1,'significant')}; % define x-tick labels
yTickData=data.idx_redcell;
y_label_text = {'.'};
x_label_text = {'Time from','stimulus (s)'};
cell_idx = 100; trial_idx = 100; % indices of example cell (plotting across all trials) or example trial (plotting across all cells)

stim.redCells= stim.unsorted(data.idx_redcell,:,:);
stim.notRedCells=stim.unsorted(data.idx_notredcell,:,:);
stim.mean.redCells= mean(stim.redCells,3); stim.mean.notRedCells = mean(stim.notRedCells,3);
stim.mean.colorGrouped= [stim.mean.redCells;stim.mean.notRedCells]; % calculates the mean across all trials for each ROI

stim.mean.normalized.colorGrouped=twoP_normMean(stim.mean.colorGrouped); % Normalize the mean to each ROI's max and min
stim.mean.normalized.sortedActivity.redCells=twoP_sortROIbyActivity(stim.mean.redCells,[]);
stim.mean.normalized.sortedActivity.notRedCells=twoP_sortROIbyActivity(stim.mean.notRedCells,[]);
stim.mean.normalized.sortedActivity.allCells=[stim.mean.normalized.sortedActivity.redCells;stim.mean.normalized.sortedActivity.notRedCells];

stim.variance.normalized.redCells = var(stim.redCells,1,3)./mean(stim.redCells,3);
stim.variance_sorted.normalized.redCells = twoP_sortROIbyActivity(stim.variance.normalized.redCells);
stim.variance.normalized.notRedCells = var(stim.notRedCells,1,3)./mean(stim.notRedCells,3);
stim.variance_sorted.normalized.notRedCells = twoP_sortROIbyActivity(stim.variance.normalized.notRedCells);
stim.variance.normalized.allCells  = [stim.variance.normalized.redCells;stim.variance.normalized.notRedCells];
stim.variance_sorted.normalized.allCells  = [stim.variance_sorted.normalized.redCells;stim.variance_sorted.normalized.notRedCells];

% singleTrial = neuralData.colorGrouped(:,:,trial_idx); % Extracts the traces of a single example trial for all neurons
%
% [sTmax,sTmaxIdx] = max(smoothdata(singleTrial,2,'sgolay',20),[],2); sTmaxIdx = sortrows([sTmaxIdx (1:1:length(sTmaxIdx))']);
% sTsorted = singleTrial(sTmaxIdx(:,2),:); % Sorts the response of all neurons during a single trial using the maximum value of the deconvolved spike trace

% Computes the mean of the responses to the side of the stimulus
[stim.mean.normalized.sideOne.redCells,stim.mean.normalized.sideTwo.redCells,stim.mean.normalized.sideOne.redCellsSideTwo,stim.mean.normalized.sideTwo.redCellsSideOne]= twoP_sortROIbyActivity(twoP_normMean(mean(stim.redCells(:,:,sideOneIdx),3)),twoP_normMean(mean(stim.redCells(:,:,sideTwoIdx),3)));
[stim.mean.normalized.sideOne.notRedCells,stim.mean.normalized.sideTwo.notRedCells,stim.mean.normalized.sideOne.notRedCellsSideTwo,stim.mean.normalized.sideTwo.notRedCellsSideOne]= twoP_sortROIbyActivity(twoP_normMean(mean(stim.notRedCells(:,:,sideOneIdx),3)),twoP_normMean(mean(stim.notRedCells(:,:,sideTwoIdx),3)));
stim.mean.normalized.sideOne.allCells = [stim.mean.normalized.sideOne.redCells;stim.mean.normalized.sideOne.notRedCells];
stim.mean.normalized.sideTwo.allCells = [stim.mean.normalized.sideTwo.redCells;stim.mean.normalized.sideTwo.notRedCells];
stim.mean.normalized.sideOne.allCellsSideTwo = [stim.mean.normalized.sideOne.redCellsSideTwo;stim.mean.normalized.sideOne.notRedCellsSideTwo];
stim.mean.normalized.sideTwo.allCellsSideOne = [stim.mean.normalized.sideTwo.redCellsSideOne;stim.mean.normalized.sideTwo.notRedCellsSideOne];

% Computes the mean of the responses to the easiest stimulus by side
[stim.mean.normalized.sideOneEasy.redCells,...
    stim.mean.normalized.sideTwoEasy.redCells,...
    stim.mean.normalized.sideOneEasy.redCellsSideTwo,...
    stim.mean.normalized.sideTwoEasy.redCellsSideOne]= ...
    twoP_sortROIbyActivity(twoP_normMean(mean(stim.redCells(:,:,trialID.leftEasyImgIdx),3)),twoP_normMean(mean(stim.redCells(:,:,trialID.rightEasyImgIdx),3)));

[stim.mean.normalized.sideOneEasy.notRedCells,...
    stim.mean.normalized.sideTwoEasy.notRedCells,...
    stim.mean.normalized.sideOneEasy.notRedCellsSideTwo,...
    stim.mean.normalized.sideTwoEasy.notRedCellsSideOne] = ...
    twoP_sortROIbyActivity(twoP_normMean(mean(stim.notRedCells(:,:,trialID.leftEasyImgIdx),3)),twoP_normMean(mean(stim.notRedCells(:,:,trialID.rightEasyImgIdx),3)));

stim.mean.normalized.sideOneEasy.allCells = [stim.mean.normalized.sideOneEasy.redCells;stim.mean.normalized.sideOneEasy.notRedCells];
stim.mean.normalized.sideTwoEasy.allCells = [stim.mean.normalized.sideTwoEasy.redCells;stim.mean.normalized.sideTwoEasy.notRedCells];
stim.mean.normalized.sideOneEasy.allCellsSideTwo = [stim.mean.normalized.sideOneEasy.redCellsSideTwo;stim.mean.normalized.sideOneEasy.notRedCellsSideTwo];
stim.mean.normalized.sideTwoEasy.allCellsSideOne = [stim.mean.normalized.sideTwoEasy.redCellsSideOne;stim.mean.normalized.sideTwoEasy.notRedCellsSideOne];

% Plot stimulus-aligned responses
hStimAlignedTraces = figure('Name',[animal ': ' session]);
set(hStimAlignedTraces, 'Units','inches',...
    'Position',[1 1 14 4]);

im_stim_aligned={stim.mean.normalized.sortedActivity.allCells,...
    stim.variance_sorted.normalized.allCells,...
    stim.mean.normalized.sideOne.allCells,...
    stim.mean.normalized.sideTwo.allCellsSideOne,...
    stim.mean.normalized.sideTwo.allCells,...
    stim.mean.normalized.sideOne.allCellsSideTwo,...
    stim.mean.normalized.sideOneEasy.allCells,...
    stim.mean.normalized.sideTwoEasy.allCellsSideOne,...
    stim.mean.normalized.sideTwoEasy.allCells,...
    stim.mean.normalized.sideOneEasy.allCellsSideTwo,...
    stim.mean.normalized.sideOneEasy.allCells-stim.mean.normalized.sideTwoEasy.allCellsSideOne,...
    stim.mean.normalized.sideTwoEasy.allCells-stim.mean.normalized.sideOneEasy.allCellsSideTwo};

imTitle={{'Mean activity','across all trials','(Grouped by cell type)'},...
    {'\sigma^{2}/\mu','across all trials'},...
    {'L stim'},...
    {'R stim','L-sorting'},...
    {'R stim',''},...
    {'L stim','R-sorting'},...
    {'L stim (DR=1)'},...
    {'R stim (DR=1)','L-sorting'},...
    {'R stim (DR=1)'},...
    {'L stim (DR=1)','R-sorting'},...
    {'L stim - R stim (DR=1)','L-sorting'},...
    {'R stim - L stim (DR=1)','R-sorting'}};
imYLabel={{'ROI #'}};
numPlotRows=4;

for i = 1:length(im_stim_aligned)
    subplot(numPlotRows,length(im_stim_aligned),[i i+length(im_stim_aligned) i+length(im_stim_aligned)*2]);
    imagesc(im_stim_aligned{i},prctile(im_stim_aligned{i}(:),[1 99])); colormap('default'); hold on;
    colormap(jet)
    hLineVert=line([find(stimWinIdx==0) find(stimWinIdx==0)], [0 size(data.neural,1)]); % Draw a line at TIME ZERO
    hLineHor=line([0 size(stim.unsorted,2)], [length(data.idx_redcell) length(data.idx_redcell)]); % Draw a line at to separate rate RED from NON-RED neurons
    set(hLineVert,'LineWidth',1.5,'Color',[1 1 1],'LineStyle',':');
    set(hLineHor,'LineWidth',1.5,'Color',[1 1 1],'LineStyle',':');
    try hTitle = title(imTitle{i}); set(hTitle,'FontName','Arial','FontSize',8,'FontWeight','bold'); end
    try
        hYLabelText=ylabel(imYLabel{i});
    end
    
    subplot(numPlotRows,length(im_stim_aligned),[i+length(im_stim_aligned)*3]);
    plot(mean(im_stim_aligned{i}));
    
    hXLabel=xlabel(x_label_text);
    set(hXLabel,'FontName','Arial','FontSize',8);
    set(gca,'YTick',yTickData ,...
        'YTickLabel', y_label_text,...
        'YColor','r',...
        'XTick',xTickStim,'XTickLabel',xTickStimLabel,'TickLength',[0 0]);
    if i>1; set(gca,'YTickLabel',[]); end
    try hYLabelText.Color=[0 0 0]; end
    
end
% print([analysisDir filesep analysisFileName '_stimulusAlignedTraces'],'-depsc');
exportgraphics(gcf,[sessAvgDir 'PSTH_' baseFileName '.pdf']);

%% Choice(lick)-aligned (all choices, including error) traces averaged across stimulus side (run previous cell first)
lickWinIdx=data.lickWinIdx; lickWinMs=data.lickWinMs; dataLick = data.dataLick;
xTickChoice=[1 find(lickWinIdx == 0) find(lickWinIdx == round(1000/data.msPerFrame)) length(lickWinIdx)]; % these tick marks are plotted
% xTickChoiceLabel = {round(lickWinMs(1)/1000,2,'significant') ,lickWinMs(lickWinMs==0), round(lickWinMs(lickWinIdx == round(1000/data.msPerFrame))/1000,2,'significant'),round(lickWinMs(end)/1000,2,'significant')};
xTickChoiceLabel = {round(lickWinMs(1)/1000,2,'significant') ,lickWinMs(lickWinMs==0), round(lickWinMs(lickWinIdx == round(1000/data.msPerFrame))/1000,2,'significant'),[]};

x_label_text = {'Time from','first lick (s)'};

% Compute the mean inferred spiking activity across all trials
meanChoice=mean(dataLick,3); % computes the mean across all trials (third dimension) for each neuron (rows); each column is a frame
nmeanChoice = twoP_normMean(meanChoice); % normalizes mean inferred spike traces for individual ROIs
mC.Sorted = twoP_sortROIbyActivity(nmeanChoice); % sorted normalized mean inferred spike traces of individual neurons

% Cross-validation analysis
mC.Odd = mean(dataLick(:,:,1:2:end),3); mC.Even = mean(dataLick(:,:,2:2:end),3);
mC.firstHalf = mean(dataLick(:,:,1:round(size(dataLick,3)/2)),3); mC.secondHalf = mean(dataLick(:,:,round(size(dataLick,3)/2)+1:end),3);
mC.norm.Odd=twoP_normMean(mC.Odd); mC.norm.Even=twoP_normMean(mC.Even);
[mC.norm.OddMax,mC.norm.OddMaxIdx]= max(mC.norm.Odd,[],2); mC.norm.OddMaxIdx= sortrows([mC.norm.OddMaxIdx (1:1:length(mC.norm.OddMaxIdx))']);
mC.norm.OddSorted = mC.norm.Odd(mC.norm.OddMaxIdx(:,2),:);
mC.norm.EvenSorted = mC.norm.Even(mC.norm.OddMaxIdx(:,2),:);

leftSideIdx = find(lick(1,lick(6,:)==1)==1); rightSideIdx = find(lick(1,lick(6,:)==1)==2);
leftSideData = dataLick(:,:,leftSideIdx); rightSideData= dataLick(:,:,rightSideIdx);
meanL=mean(leftSideData,3); meanR=mean(rightSideData,3);
nmeanL=twoP_normMean(meanL); nmeanR=twoP_normMean(meanR);
[mC.norm.LL,mC.norm.RR,mC.norm.LR,mC.norm.RL]=twoP_sortROIbyActivity(nmeanL,nmeanR);

[easyLeftTrialID,iB,iEasyLeft] = intersect(find(ratioEvents==0)',data.dataLickTrialNumbers); [easyRightTrialID,iB,iEasyRight] = intersect(find(ratioEvents==1)',data.dataLickTrialNumbers); 
dataEasyLeft=dataLick(:,:,iEasyLeft); dataEasyRight=dataLick(:,:,iEasyRight);
meanEasyL=mean(dataEasyLeft,3); meanEasyR=mean(dataEasyRight,3);
nmeanLeasy=twoP_normMean(meanEasyL); nmeanReasy=twoP_normMean(meanEasyR);
[mC.norm.LLeasy,mC.norm.RReasy,mC.norm.LReasy,mC.norm.RLeasy]=twoP_sortROIbyActivity(nmeanLeasy,nmeanReasy);

hStimSideTraces = figure('Name',[animal ': ' session]);
set(hStimSideTraces, 'Units','inches',...
    'Position',[1 1 10 4]);

im={mC.Sorted, mC.norm.OddSorted, mC.norm.EvenSorted, mC.norm.LL, mC.norm.RL, mC.norm.RR, mC.norm.LR, mC.norm.LLeasy, mC.norm.RLeasy, mC.norm.RReasy, mC.norm.LReasy};
imTitle={{'Choice-aligned activity','all correct trials'},...
    {'Correct choices -','odd trials'},...
    {'Correct choices -','even trials','(odd-sorted)'},...
    {'Left choices', 'left sorted'},...
    {'Right choices','left sorted'},...
    {'Right choices','right sorted)'},...
    {'Left choices','right sorted'},...
    {'Left choices - easy', 'left sorted'},...
    {'Right choices - easy','left sorted'},...
    {'Right choices - easy','right sorted)'},...
    {'Left choices - easy','right sorted'}};
imYLabel={{'ROI #'}};
numPlotRows=1;

for i = 1:numel(im)
    subplot(numPlotRows,numel(im),i);
    imagesc(im{i},prctile(im{i}(:),[10 90])); colormap('default'); hold on;
    hLineVert=line([find(lickWinIdx==0) find(lickWinIdx==0)], [0 size(dataLick,1)]); % Draw a line at TIME ZERO
    set(hLineVert,'LineWidth',1.5,'Color',[1 1 1],'LineStyle',':');
    try hTitle = title(imTitle{i}); set(hTitle,'FontName','Arial','FontSize',8,'FontWeight','bold'); end
    try
        ylabel(imYLabel{i});
    end
    
    hXLabel=xlabel(x_label_text);
    set(hXLabel,'FontName','Arial','FontSize',8);
    set(gca,'XTick',xTickChoice,'XTickLabel',xTickChoiceLabel,'TickLength',[0 0]);
    if i>1; set(gca,'YTickLabel',[]); end
end

% print([analysisDir filesep analysisFileName '_choiceAlignedTraces'],'-depsc');
exportgraphics(gcf,[sessAvgDir 'PETH_Lick_' baseFileName '.pdf']);

%% 
iStart=-3; iEnd = -1; % index (relative to choice onset at index 0) to be averaged into a time bin for analysis
X = squeeze(mean(dataLick(:,find(lickWinIdx==iStart):find(lickWinIdx==iEnd),:),2))';
Y = lick(1,:)';

% permROI = randperm(size(dataLick,1));
% Xone = squeeze(mean(dataLick(permROI(1:floor(length(permROI)/2)),find(lickWinIdx==-3):find(lickWinIdx==-1),:),2))';
% Xtwo = squeeze(mean(dataLick(permROI(floor(length(permROI)/2)+1:end),find(lickWinIdx==-3):find(lickWinIdx==-1),:),2))';
for i = 1:length(lick(1,:))
    if lick(1,i) == 1
        Y{i,1}='leftChoice';
    elseif lick(1,i) == 2
        Y{i,1} = 'rightChoice';
    end
end

% SVMModel = fitclinear(Xone,Y,'ClassNames',{'leftChoice','rightChoice'},...
%     'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',struct('MaxObjectiveEvaluations',50,'AcquisitionFunctionName','expected-improvement-plus','KFold',40,'UseParallel',true));

SVMModel = fitcsvm(X,Y,'ClassNames',{'leftChoice','rightChoice'},...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus','KFold',10,'UseParallel',false ));


% Ypred_self=predict(SVMModel,Xone); 
% Ypred_other=predict(SVMModel,Xtwo); 

% disp(['Accuracy in predicting self: ' num2str((sum(strcmp(Y,Ypred_self))/numel(strcmp(Y,Ypred_self))*100)) '%'])
% disp(['Accuracy in predicting other: ' num2str((sum(strcmp(Y,Ypred_other))/numel(strcmp(Y,Ypred_other))*100)) '%'])

% numIters=1000;
% Ypred_all = zeros(1,numIters);
% parfor i = 1:numIters
%     permROI = randperm(size(dataLick,1));
%     Xtwo = squeeze(mean(dataLick(permROI(floor(length(permROI)/2)+1:end),find(lickWinIdx==-3):find(lickWinIdx==-1),:),2))';
%     Ypred_two=predict(SVMModel,Xtwo); 
%     Ypred_all(i)=sum(strcmp(Y,Ypred_two))/numel(strcmp(Y,Ypred_two));
% end
% histogram(Ypred_all,100);

% SVMModel = fitrlinear(X,Y,...
%     'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',struct('MaxObjectiveEvaluations',1200,'AcquisitionFunctionName','expected-improvement-plus','KFold',40,'UseParallel',true));

% n1=1; n2=2;
% n1L= mean(squeeze(dataLick(n1,find(lickWinIdx==-6):find(lickWinIdx==-1),lick(1,:)==1)));
% n1R= mean(squeeze(dataLick(n1,find(lickWinIdx==-6):find(lickWinIdx==-1),lick(1,:)==2)));
% n2L= mean(squeeze(dataLick(n2,find(lickWinIdx==-6):find(lickWinIdx==-1),lick(1,:)==1)));
% n2R= mean(squeeze(dataLick(n2,find(lickWinIdx==-6):find(lickWinIdx==-1),lick(1,:)==2)));
% 
% figure(3);
% plot(n1L,n2L,'.b'); hold on; plot(n1R,n2R,'.r')