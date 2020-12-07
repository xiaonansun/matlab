% function twoP_scripts(animal,session)
% mfilename('fullpath')
%% Load trial-aligned fluorescence data
clear all; close all;

% Specify session to load
animal = 'Fez51';
session = '20200402a';
analysisFileName = [animal '_' session];
addpath(genpath('C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab\2pAnalysisCode'));

% Define directories
if convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('MANHASSET')
    imagingRootDir = 'G:\2PData';
else
%     imagingRootDir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
    imagingRootDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
%     imagingRootDir = 'F:\suite2p';
    disp('Loading suite2p data...');
end
imagingSubDir = 'suite2p\plane0';

disp(['Animal: ' animal '; Session: ' session]);

% --- Load the google sheets document "2photon acquisition record" --- %
docid = '1ADcwZJygK7fV0zOq537V1Vsyv45-OF3WkMopNbjFifY';
expTable=GetGoogleSpreadsheet(docid); % this function (GetGoogleSpreadsheet.m) needs to be downloaded
bhvColIdx=find(contains(expTable(1,:),'Behavior file name'));
iFolderColIdx=find(contains(expTable(1,:),'Folder'));
try
bhvRowIdx=find(contains(expTable(:,iFolderColIdx),session));
bhvFName=expTable{bhvRowIdx(contains(expTable(bhvRowIdx),animal)),bhvColIdx};
catch ME
    disp([ME.message]);
    disp('Cannot load session. Please check the session name input.');
    analysis.error.behaviorTable = ME.message;
end

analysisDir = 'C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab\twoP_analysis'; %directory of the analysis history summary file
load([analysisDir filesep 'analysis.mat']); % load the analysis history summary file
currentDate=datestr(now,'yyyy-mm-dd_HHMMss'); currAnalIdx=length(analysis)+1; % append data of the current analysis to the summary file

fprintf('Loading 2P imaging data...');
npy = twoP_importSuite2pData(animal,session, imagingRootDir); % organize suite2p data into npy files (ROI and spiking)
try
data = twoP_alignDetectionTask(npy.ops, npy.spks, npy.iscell, npy.redcell, npy.bin_MScan_filepath); % align suite2p data to sensory stimulus
catch ME
end

if exist('ME','var')
    disp(['Error detected: ' ME.identifier '. ' ME.message])
else
    disp('2P DATA LOADED!');
end

analysis(currAnalIdx).date=currentDate;
analysis(currAnalIdx).animal=animal; analysis(currAnalIdx).session=session;
save([analysisDir filesep 'analysis.mat'],'analysis'); % save a historical record of this analysis script that has been run
writetable(struct2table(analysis),[analysisDir filesep 'analysis.csv'])

% Load behavior data
fprintf('Loading behavior data...');
[SessionData,bhvFilePath] = twoP_loadBehaviorSession(animal,session,bhvFName); SessionData = SessionData.SessionData;
% [filepath,bhvFileName,ext]=fileparts(bhvFilePath); clear filepath ext;
fprintf('DONE!\n');
% allStim_R_L = SessionData.CorrectSide;

disp(['Directory of suite2p output: ' fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir)]);
disp(['Path of behavior data: ' bhvFilePath]);
disp(['Number of trials in neural data matrix: ' num2str(length(data.trialNumbers))]);
disp(['Number of trials recorded by 2P imaging: ' num2str(max(data.trialNumbers))]);
disp(['Number of trials recorded by Bpod: ' num2str(length(SessionData.Rewarded))]);
if abs(max(data.trialNumbers) -length(SessionData.Rewarded)) >= 10; disp('Behavior and imaging differs by more than 10 trials: do you have the correct behavior file?'); end;

% save(fullfile(imagingRootDir,animal,session,imagingSubDir,'dataRocAnalysis.mat'),'firstSideTryAl','allStim_R_L','SessionData');
% B = regexp(s1,'\d*','Match');
%% makes adjustments to data struct

% Discard extra trials recorded by the two-photon, usually this is one extra trial at the end of the session
if length(data.trialNumbers) < size(data.neural,3) 
    disp(['Note: there are ' num2str(size(data.neural,3) - length(data.trialNumbers)) ' fewer trial ID(s) compared to the number of trials in the neural data array.'])
    delIdx = find(~ismember(1:1:size(data.neural,3), 1:1:length(data.trialNumbers)));
    data.neural(:,:,delIdx)=[];
    data.stimFramesOrig(delIdx)=[];
    data.DS(:,:,delIdx)=[];
    data.stimSamplesOrig(delIdx)=[];
    data.analog(:,:,delIdx)=[];
end
% Compares the number of imaging trials to the the number of trials recorded in the behavior data and discard extra imaging trials
if max(data.trialNumbers) > length(SessionData.Rewarded) 
    disp(['Note: ' num2str(max(data.trialNumbers) - length(SessionData.Rewarded)) ' additional trial(s) were recorded during imaging than during behavior training.'])
    data.trialNumbersOrig = data.trialNumbers;
    delIdx = find(~ismember(data.trialNumbers, 1:1:length(SessionData.Rewarded)));
    data.trialNumbers(delIdx)=[];
    data.neural(:,:,delIdx)=[];
    data.stimFramesOrig(delIdx)=[];
    data.DS(:,:,delIdx)=[];
    data.stimSamplesOrig(delIdx)=[];
    data.analog(:,:,delIdx)=[];
end

%% Define event-aligned matrices
fps=1000/data.msPerFrame;
stdThresh=2.5; % Reject delayed lick responses above a threshold, expressed in units of standard deviations from the mean
filterWin = 600; % length of smoothing filter in milliseconds
filterFrames = round(filterWin/data.msPerFrame); % number of frames used for the smoothing filter as specified by the time constant in line above
iTrialID=data.trialNumbers;
[stimTime, spoutTime, lickR, lickL, water]=behavior_findEventTiming(SessionData); %Use the findEventTiming function to extract the timing of stimulus onset, spout movement, and licks.
stimTime=stimTime(data.trialNumbers); spoutTime=spoutTime(data.trialNumbers);
if SessionData.nTrials > length(water); water=[water zeros(1,SessionData.nTrials-length(water))]; end
water=water(data.trialNumbers);
if SessionData.nTrials > length(lickR); lickR = [lickR num2cell(zeros(1,SessionData.nTrials-length(lickR)))]; end
if SessionData.nTrials > length(lickL); lickL = [lickL num2cell(zeros(1,SessionData.nTrials-length(lickL)))]; end

lickR=lickR(data.trialNumbers); lickL=lickL(data.trialNumbers);

stimOnsetIdx=find(data.neuralTimes==0);
minLicks=2;

lick= zeros(3,length(data.trialNumbers));

% Generate a row of response side: 1 is left and 2 is right
lick(1,:) = SessionData.ResponseSide(data.trialNumbers);

% Generate a row of the timing (in seconds) of the first lick event
% relative to trial onset
leftLickIdx=find(cellfun(@numel,lickL)>=minLicks);
rightLickIdx=find(cellfun(@numel,lickR)>=minLicks);
[val, ir, il]=intersect(rightLickIdx,leftLickIdx);
rightLickIdx(ir)=[]; leftLickIdx(il)=[]; % remove trials where the animal licked both spouts more than once
leftLickTiming=cellfun(@(x)x(1),lickL((leftLickIdx)));
rightLickTiming=cellfun(@(x)x(1),lickR(rightLickIdx)); %compute the timing of left or right licks
lick(2,leftLickIdx)=leftLickTiming; lick(2,rightLickIdx)=rightLickTiming; % compute the timing of the first lick of a true lick response
lick(2,lick(2,:)==0)=NaN;
lick(2,lick(2,:)>nanmean(lick(2,:))+nanstd(lick(2,:))*stdThresh)=NaN;

% compute the timing of the first lick of a true lick response relative to
% the onset of the sensory stimulus
lick(3,:)=lick(2,:)-stimTime;

% Compute the time index of the first lick
lick(4,:)= round(lick(3,:)*fps)+stimOnsetIdx;
maxPostLickIdx=max(length(data.neuralTimes)-lick(4,:));

lick(5,:)= data.trialNumbers; % Trial ID
lick(6,:)=SessionData.Rewarded(lick(5,:));
rewardedTrials = lick(5,lick(6,:)==1);
lick(:,isnan(lick(1,:)))=[];
[row,col,v] = find(isnan(lick)); lick(:,col)=[];

preLickFrames = round(fps); postLickFrames = size(data.neural,2)-max(lick(4,:)); % 1 second windows before stimulus
lickWinIdx = -preLickFrames:1:postLickFrames; % IMPORTANT: frame index relative to animal lick
lickWinMs = lickWinIdx*data.msPerFrame;  % IMPORTANT: time (in milliseconds) relative to animal lick
startIdx = lick(4,:)'- preLickFrames -1 ;  endIdx = lick(4,:)' + postLickFrames;
dataLick = zeros(size(data.neural,1),endIdx(1)-startIdx(1),length(lick));
[C, iL, iD]=intersect(lick(5,:),data.trialNumbers);
for i = 1:length(lick); dataLick(:,:,i) = data.neural(:,startIdx(i):endIdx(i)-1,iD(i)); end

% histogram(rightTiming, 50); hold on; histogram(leftTiming,50);
%% Behavior analysis
events =SessionData.stimEvents;
% organize stimulus events of each trial (row) into columns: column 1 is
% left and column 2 is right
events = cellfun(@(x) x{:},[cellfun(@(x) x(1),events,'UniformOutput',false)' cellfun(@(x) x(2),events,'UniformOutput',false)'],'UniformOutput',false);
frEvents= cellfun(@numel,events);
frEvents(:,3)= frEvents(:,1)./sum(frEvents,2); % ratio of left to all events presented during a trial; 0=0% left-sided events, 1=100% left-sided events
frEvents(:,4)=(SessionData.CorrectSide)';
frEvents(:,5)=(SessionData.DistStim)'; % distractor

%% Reorganize data for PCA and dPCA
% trialID.imaging is the master trial index of 2p imaging data. This
% structure field should contain all indices contained by the other trialID
% fields defined here

data.neuralNorm=twoP_normMean(data.neural);
time=data.neuralTimesMs;
numOfBins = 6;
edges = 0:1/numOfBins:1; % divide the auditory stimuli presented throughout the session into bins which in this case is 11: 0
edges = [0 0.01 edges(2:end-1) 0.99 1];

stimValues=frEvents(:,3); % fraction of left sided events

% Remove auto-reward trials
trialID.notAutoReward=find(~SessionData.AutoReward); % Find the IDs of non-auto-reward trials
[trialID.imaging,trialID.trialNumbersIdx,trialID.AutoRewardIdx]= intersect(data.trialNumbers,trialID.notAutoReward); % Find the non-auto-reward imaging trial IDs

% ADDITIONALLY, WILL ALSO NEED TO REMOVE TRAINING MODE TRIALS: SessionData.SingleSpout

% remove
trialID.leftLick = find(SessionData.ResponseSide==1); trialID.rightLick = find(SessionData.ResponseSide==2); % Find the behavior trial IDs of left or right decisions
trialID.imagingLeftLick = intersect(trialID.imaging,trialID.leftLick); % find the trial IDs of non-auto-reward left decisions
trialID.imagingRightLick = intersect(trialID.imaging,trialID.rightLick); % find the trial IDs of non-auto-reward right decisions

[binCounts,edges,stimulusCondition]=histcounts(stimValues(trialID.imaging),edges); % divides the
trialID.stimulusCondition = cell(max(unique(stimulusCondition)),1); % creates a field of the cell class to hold trial ID information for each stimulus category

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

% dimensions: (1) neuron ID, (2) time, (3) trials grouped by stimulus
% condition, (4) stimulus condition, (5) decision: left or right
data.pcaArray=nan(size(data.neuralNorm,1),size(data.neuralNorm,2),maxStimTrials,length(trialID.stimulusCondition),2);
m=0;
for i = 1:size(data.pcaArray,5) % i: response (left or right)
    for j = 1:size(data.pcaArray,4) % j: stimulus category (1 through 11 with 1=left? and 11=being right? and everything else in between)
        for k = 1:length(stimConditions{i}{j})
            data.pcaArray(:,:,k,j,i) = data.neuralNorm(:,:,ismember(data.trialNumbers,stimConditions{i}{j}(k)));
            m=m+1;
        end
    end
end

% Converts the array to new dimensions: (1) neuron ID, (2) stimulus conditions, (3) Decision (L/R) ouped by stimulus
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

%% Re-organize data for plotting tuning curves
data.neuralNorm=twoP_normMean(data.neural);
time=data.neuralTimesMs;
stimRange = [500 1000]; % stimulus start and end time (in milliseconds)
idxtime=find(time > stimRange(1) & time <= stimRange(2));
numOfBins = 5; % Number of bins (to be specified by the user) excluding the two extremes: 0 and 1. In total, the number of bins will be numOfBins + 2
edges = 0:1/numOfBins:1; % divide the auditory stimuli presented throughout the session into bins 
edges = [0 0.01 edges(2:end-1) 0.99 1];

stimValues=frEvents(:,3); % fraction of left sided events
[binCounts,edges,stimulusCondition]=histcounts(stimValues(trialID.imaging),edges); % divides the vector (each entry is a trial) containing the fraction of left-sided events into bins
sR=nan(size(data.neural,1),length(idxtime),max(binCounts),length(binCounts)); % sR = stimulus response. Creates a 4-D vector. 1: cell ID, 2: time, 3: trial, 4: stimulus condition (from all left- to all right-sided events)

for i = 1:length(binCounts)
% sR(:,:,1:binCounts(i),i) = data.neural(:,idxtime,trialID.imaging(stimulusCondition==i));
sR(:,:,1:binCounts(i),i) = data.neural(:,idxtime,stimulusCondition==i);
end

tM = squeeze(nanmean(sum(sR,2),3)); % tM = tuning matrix, where rows are cells and columns are the stimulus conditions from 100% left-sided events to 100% right-sided events.
[sMax,sMaxIdx]=max(tM,[],2);
sMaxIdx(:,2)=1:length(sMaxIdx); sMaxIdx(:,3:4)=sortrows(sMaxIdx(:,1:2),1); %sort the matrix by the stimulus condition with the maximum response
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



%% Plot stimulus-aligned mean (across trials) traces
filterType = 'sgolay';
stimWinIdx= data.neuralTimes;
trialID.allImaging=lick(5,:);
trialID.all=intersect(trialID.allImaging,find(~SessionData.AutoReward));
trialID.leftEasy = intersect(trialID.all,find(frEvents(:,3)==1));
[trialID.leftEasyIdx, trialID.leftEasySubIdx,trialID.leftEasyImgIdx]= intersect(trialID.leftEasy,trialID.allImaging);
trialID.rightEasy = intersect(trialID.all,find(frEvents(:,3)==0));
[trialID.rightEasyIdx,trialID.rightEasySubIdx,trialID.rightEasyImgIdx] = intersect(trialID.rightEasy,trialID.allImaging);
correctSide = SessionData.CorrectSide(trialID.all);
sideOneIdx = find(correctSide==1); sideTwoIdx = find(correctSide==2);
stim.unsorted= data.neural(:,:,iD);

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
stim.mean.normalized.sortedActivity.redCells=twoP_sortROIbyActivity(stim.mean.redCells,[],filterType,filterWin/fps);
stim.mean.normalized.sortedActivity.notRedCells=twoP_sortROIbyActivity(stim.mean.notRedCells,[],filterType,filterWin/fps);
stim.mean.normalized.sortedActivity.allCells=[stim.mean.normalized.sortedActivity.redCells;stim.mean.normalized.sortedActivity.notRedCells];

stim.variance.normalized.redCells = var(stim.redCells,1,3)./mean(stim.redCells,3);
stim.variance.normalized.notRedCells = var(stim.notRedCells,1,3)./mean(stim.notRedCells,3);
stim.variance.normalized.allCells  = [stim.variance.normalized.redCells;stim.variance.normalized.notRedCells];

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
hStimAlignedTraces = figure(1);
set(hStimAlignedTraces, 'Units','inches',...
    'Position',[1 1 14 4]);

im_stim_aligned={stim.mean.normalized.sortedActivity.allCells,...
    stim.variance.normalized.allCells,...
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
    hLineVert=line([find(stimWinIdx==0) find(stimWinIdx==0)], [0 size(dataLick,1)]); % Draw a line at TIME ZERO
    hLineHor=line([0 size(stim.unsorted,2)], [length(data.idx_redcell) length(data.idx_redcell)]); % Draw a line at TIME ZERO
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

print([analysisDir filesep analysisFileName '_stimulusAlignedTraces'],'-depsc');

%% Choice(lick)-aligned (all choices, including error) traces averaged across stimulus side (run previous cell first)
xTickChoice=[1 find(lickWinIdx == 0) find(lickWinIdx == round(1000/data.msPerFrame)) length(lickWinIdx)]; % these tick marks are plotted
xTickChoiceLabel = {round(lickWinMs(1)/1000,2,'significant') ,lickWinMs(lickWinMs==0), round(lickWinMs(lickWinIdx == round(1000/data.msPerFrame))/1000,2,'significant'),round(lickWinMs(end)/1000,1,'significant')};
x_label_text = {'Time from','first lick (s)'};

% Computing the mean inferred spiking activity across all trials
meanChoice=mean(dataLick,3); % computes the mean across all trials (third dimension) for each neuron (rows); each column is a frame
nmeanChoice=bsxfun(@rdivide,bsxfun(@minus,meanChoice,min(meanChoice,[],2)),max(meanChoice,[],2)-min(meanChoice,[],2)); % normalizes the mean inferred (deconvolved) spike trace
[NMmaxChoice,NMmaxChoiceIdx]= max(smoothdata(nmeanChoice,2,'sgolay',filterFrames),[],2); NMmaxChoiceIdx= sortrows([NMmaxChoiceIdx (1:1:length(NMmaxChoiceIdx))']); % (1) smoothes the normalized mean inferred spike trace and (2) sorts neurons based on when the maximum spike rate occurs
mC.Sorted = nmeanChoice(NMmaxChoiceIdx(:,2),:); % this is the sorted normalized mean inferred spike traces of individual neurons

% Cross-valication analysis
mC.Odd = mean(dataLick(:,:,1:2:end),3); mC.Even = mean(dataLick(:,:,2:2:end),3);
mC.firstHalf = mean(dataLick(:,:,1:round(size(dataLick,3)/2)),3); mC.secondHalf = mean(dataLick(:,:,round(size(dataLick,3)/2)+1:end),3);
mC.norm.Odd=bsxfun(@rdivide,bsxfun(@minus,mC.Odd,min(mC.Odd,[],2)),max(mC.Odd,[],2)-min(mC.Odd,[],2));
mC.norm.Even=bsxfun(@rdivide,bsxfun(@minus,mC.Even,min(mC.Even,[],2)),max(mC.Even,[],2)-min(mC.Even,[],2));
[mC.norm.OddMax,mC.norm.OddMaxIdx]= max(smoothdata(mC.norm.Odd,2,'sgolay',filterFrames),[],2); mC.norm.OddMaxIdx= sortrows([mC.norm.OddMaxIdx (1:1:length(mC.norm.OddMaxIdx))']);
mC.norm.OddSorted = mC.norm.Odd(mC.norm.OddMaxIdx(:,2),:);
mC.norm.EvenSorted = mC.norm.Even(mC.norm.OddMaxIdx(:,2),:);

leftSideIdx = find(lick(1,lick(6,:)==1)==1); rightSideIdx = find(lick(1,lick(6,:)==1)==2);
leftSideData = dataLick(:,:,leftSideIdx); rightSideData= dataLick(:,:,rightSideIdx);
meanL=mean(leftSideData,3); meanR=mean(rightSideData,3);
nmeanL=bsxfun(@rdivide,bsxfun(@minus,meanL,min(meanL,[],2)),max(meanL,[],2)-min(meanL,[],2));
nmeanR=bsxfun(@rdivide,bsxfun(@minus,meanR,min(meanR,[],2)),max(meanR,[],2)-min(meanR,[],2));
[NMmaxL,NMmaxLIdx]= max(smoothdata(nmeanL,2,'sgolay',filterFrames),[],2); NMmaxLIdx= sortrows([NMmaxLIdx (1:1:length(NMmaxLIdx))']);
[NMmaxR,NMmaxRIdx]= max(smoothdata(nmeanR,2,'sgolay',filterFrames),[],2); NMmaxRIdx= sortrows([NMmaxRIdx (1:1:length(NMmaxRIdx))']);
mC.norm.LL = nmeanL(NMmaxLIdx(:,2),:); % plots LEFT choices sorted by LEFT side
mC.norm.RL = nmeanR(NMmaxLIdx(:,2),:); %plots RIGHT choices sorted by LEFT side
mC.norm.RR = nmeanR(NMmaxRIdx(:,2),:); % plots RIGHT choices sorted by RIGHT side
mC.norm.LR = nmeanL(NMmaxRIdx(:,2),:); % plots LEFT choices sorted by RIGHT side

hStimSideTraces = figure(2);
set(hStimSideTraces, 'Units','inches',...
    'Position',[1 1 10 4]);

im={mC.Sorted, mC.norm.OddSorted, mC.norm.EvenSorted, mC.norm.LL, mC.norm.RL, mC.norm.RR, mC.norm.LR};
imTitle={{'Choice-aligned activity','all correct trials'},...
    {'Correct choices -','odd trials'},...
    {'Correct choices -','even trials','(odd-sorted)'},...
    {'Left choices', 'left sorted'},...
    {'Right choices','left sorted'},...
    {'Right choices','right sorted)'},...
    {'Left choices','right sorted'}};
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

print([analysisDir filesep analysisFileName '_choiceAlignedTraces'],'-depsc');

%% Plot discrimination curve
hDisc = figure(3);
cInd = 1:length(SessionData.Rewarded);
[distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(SessionData, cInd);
% rateDisc_audioDiscCurve(bhv, cInd, distBins, discOnly, fixBias, returnCIs)
plot(distRatio,pChoseHigh,'.-k');
xlim([0 1]); ylim([0 1]);
title([num2str(sum(SessionData.Rewarded)) ' trials rewarded.']);

%%
Animal = 'Fez73';
cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
% [Performance,bhv,allDates] = RateDisc_learningCurves(Animal,cPath);
[Performance,bhv] = DelayedLoc_learningCurves(Animal,cPath);
