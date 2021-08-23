function Performance = SpatialDisc_ReverseCorrelation(Animal,path,binSize,stimDuration,lSessions,highDetection,mmOnly)
% Analyze behavioral data from SpatialDisc paradigm to get the psychophysical reverse
% correlation. Reads all available files in 'path', ignoring data files that contain less then 100 trials.
% 
% Inputs:
% Animal: Name of animal.
% path: Path to behavioral data without separator at the end. E.g. 'W:\data\Behavior_Simon'
% 
% Optional inputs: 
% binSize: Size of timebins for reverse correlation. Default is 10ms.
% stimDuration: Duration of a single trial in seconds. Default is 1.
% lSessions: Only use last 'lSessions' sessions for behavioral analysis.
% highDetection: Only use sessions at which detection was at 90% or higher. This is the default.
% mmOnly: Only use sessions at which multisensory stimuli were presented. This is the default.

%% check optional input
if ~exist('binSize','var') || isempty(binSize)
    binSize = 0.01;
else
    binSize = binSize/1000; % convert to seconds.
end

if ~exist('stimDuration','var') || isempty(stimDuration)
    stimDuration = 1;
end

if ~exist('lSessions','var') || isempty(lSessions)
    lSessions = inf;
end

if ~exist('highDetection','var') || isempty(highDetection)
    highDetection = true;
end

if ~exist('mmOnly','var') || isempty(mmOnly)
    mmOnly = true;
end

%% run basic analysis to create appended array
[Performance,bhv] = SpatialDisc_BasicAnalysis(Animal{1},path,lSessions,highDetection,mmOnly);

%% compute reverse correlation
leftRate = zeros(length(bhv.stimEvents),round(stimDuration/binSize)); %allocate leftRate array
rightRate = zeros(length(bhv.stimEvents),round(stimDuration/binSize)); %allocate rightRate array

for iTrials = 1:length(bhv.stimEvents)
    leftRate(iTrials,:) = histcounts(bhv.stimEvents{iTrials}{3},0:binSize:stimDuration);  %count events in each bin for left events
    rightRate(iTrials,:) = histcounts(bhv.stimEvents{iTrials}{4},0:binSize:stimDuration);  %count events in each bin for right events
end

excessRate = leftRate-rightRate; %compute excess rate
distFractions = bhv.DistStim ./ bhv.TargStim;
lInd = bhv.Assisted & ((bhv.CorrectSide == 1 & bhv.Rewarded) |(bhv.CorrectSide == 2 & bhv.Punished)); %index for non-detection, left-choice trials
rInd = bhv.Assisted & ((bhv.CorrectSide == 2 & bhv.Rewarded) |(bhv.CorrectSide == 1 & bhv.Punished)); %index for non-detection, right-choice trials

usedRate = 0;
trialSelect = distFractions >= usedRate;
lWeights = mnrfit(excessRate(trialSelect,ceil(0.05/binSize)+1:end-ceil(0.05/binSize)),lInd(trialSelect)'+1);
rWeights = mnrfit(excessRate(trialSelect,ceil(0.05/binSize)+1:end-ceil(0.05/binSize)),rInd(trialSelect)'+1);
lWeights = [zeros(round(0.05/binSize),1);lWeights(2:end)];lWeights(1/binSize) = 0;
rWeights = [zeros(round(0.05/binSize),1);rWeights(2:end)];rWeights(1/binSize) = 0;

figure;
plot(binSize:binSize:1,smooth(lWeights),'g','linewidth',2);hold on
plot(binSize:binSize:1,smooth(rWeights),'r','linewidth',2); 
title([Animal ' - logistic regression weights']);
axis square; hline(0); 
ylabel('beta weights');xlabel('time(s)');
% ylim([-.5 .5]);
end