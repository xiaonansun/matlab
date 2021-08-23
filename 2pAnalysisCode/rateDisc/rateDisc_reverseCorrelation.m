function rateDisc_reverseCorrelation(bhv, groupnames, binSize, stimDur)
% Analyze behavioral data from SpatialDisc paradigm to get the psychophysical reverse
% correlation.
%
% Optional inputs:
% binSize: Size of timebins for reverse correlation. Default is 10ms.
% stimDur: Duration of a single trial in seconds. Default is 1.

%% check optional input
if ~exist('binSize','var') || isempty(binSize)
    binSize = 0.125;
end

if ~exist('stimDur','var') || isempty(stimDur)
    stimDur = 1;
end

stimDur = round(stimDur/binSize) * binSize; %make sure stimDur is can be divided by binSize
fiberLocations = { 'Frontal' 'Parietal'};
fiberColors = {[0 0 1] [1 0 0]};

%%
for cLoc = 1 : 2
    figure; Cnt = 0;
    for x = 1 : length(bhv)
        %% compute reverse correlation
        oInd = (bhv{x}.optoDur == 1.3 | bhv{x}.optoDur == 1.5 | bhv{x}.optoDur == 0.6) & bhv{x}.optoSide == 3 & bhv{x}.optoPower > 1; %bilateral stim and delay stimuation at max. power.
        pInd = ~bhv{x}.DidNotChoose & ~bhv{x}.DidNotLever & logical(bhv{x}.Assisted & bhv{x}.DistStim > 0); %only use active discrimination trials
        
        cInd = pInd & oInd & bhv{x}.stimLocation == cLoc & bhv{x}.optoType == 1; %optogenetic; frontal stim only
        [lOptoWeights, rOptoWeights] = reverseCorrelation(cInd, bhv{x}, binSize, stimDur); %get optogenetic weights
        
        cInd = pInd & bhv{x}.optoDur == 0 & ismember(bhv{x}.SessionNr, unique(bhv{x}.SessionNr(cInd))); %non-optogenetic trials from same sessions
        [lCtrlWeights, rCtrlWeights] = reverseCorrelation(cInd, bhv{x}, binSize, stimDur); %get optogenetic weights
        
        %% show figure for weights
        Cnt = Cnt+1;
        subplot(length(bhv),3,Cnt)
        stdshade(lCtrlWeights',0:binSize:stimDur,'g',0.5,3); hold on
        stdshade(rCtrlWeights',0:binSize:stimDur,'r',0.5,3);
        title([groupnames{x} ' - Control weights']);
        axis square; hline(0); ylim([-0.12 0.12]);
        ylabel('beta weights');xlabel('time(s)');
        
        Cnt = Cnt+1;
        subplot(length(bhv),3,Cnt)
        stdshade((rCtrlWeights - lCtrlWeights)',0:binSize:stimDur,'k',0.5,3); hold on
        stdshade((rOptoWeights - lOptoWeights)',0:binSize:stimDur,fiberColors{cLoc},0.5,3);
        title([fiberLocations{cLoc} ' - L/R weight difference']);
        axis square; hline(0); ylim([-0.05 0.25]);
        ylabel('beta weights');xlabel('time(s)');
        
        Cnt = Cnt+1;
        subplot(length(bhv),3,Cnt)
        stdshade((rCtrlWeights - lCtrlWeights)' - (rOptoWeights - lOptoWeights)',0:binSize:stimDur,'k',0.5,3); hold on
        title([fiberLocations{cLoc} ' - control/opto difference']);
        axis square; hline(0); ylim([-0.1 0.2]);
        ylabel('beta weights');xlabel('time(s)'); drawnow;
    end
end



function [lWeights, rWeights] = reverseCorrelation(cInd, cBhv, binSize, stimDur)
%% run logistic regresion
lWeights = NaN(stimDur/binSize+1, length(unique(cBhv.AnimalID)) + 1);
rWeights = NaN(stimDur/binSize+1, length(unique(cBhv.AnimalID)) + 1);
for iAnimals = 1 : length(unique(cBhv.AnimalID)) + 1
    % stimDur = round(max(bhv.stimDur)/binSize) * binSize; %make sure stimDur is can be divided by binSize
    leftRate = zeros(length(cBhv.stimEvents),round(stimDur/binSize)); %allocate leftRate array
    rightRate = zeros(length(cBhv.stimEvents),round(stimDur/binSize)); %allocate rightRate array
    
    for iTrials = 1:length(cBhv.stimEvents)
        leftRate(iTrials,:) = histcounts(cBhv.stimEvents{iTrials}{1},0:binSize:stimDur);  %count events in each bin for left events
        rightRate(iTrials,:) = histcounts(cBhv.stimEvents{iTrials}{2},0:binSize:stimDur);  %count events in each bin for right events
    end
    
    excessRate = leftRate-rightRate; %compute excess rate
    distFractions = cBhv.DistStim ./ cBhv.TargStim;
    lInd = cBhv.Assisted & ((cBhv.CorrectSide == 1 & cBhv.Rewarded) |(cBhv.CorrectSide == 2 & cBhv.Punished)); %index for non-detection, left-choice trials
    rInd = cBhv.Assisted & ((cBhv.CorrectSide == 2 & cBhv.Rewarded) |(cBhv.CorrectSide == 1 & cBhv.Punished)); %index for non-detection, right-choice trials
    
    if iAnimals > length(unique(cBhv.AnimalID))
        trialSelect = cInd;
    else
        trialSelect = cInd & cBhv.AnimalID == iAnimals;
    end
    
    if sum(trialSelect) > 100
        lWeights(:, iAnimals) = mnrfit(excessRate(trialSelect,:),lInd(trialSelect)'+1);
        rWeights(:, iAnimals) = mnrfit(excessRate(trialSelect,:),rInd(trialSelect)'+1);
        lWeights(:, iAnimals) = [zeros(ceil(0.05/binSize),1);lWeights(ceil(0.05/binSize)+1:end,iAnimals)]; %first 50ms shouldnt be used due to deadtime after onset pulse
        rWeights(:, iAnimals) = [zeros(ceil(0.05/binSize),1);rWeights(ceil(0.05/binSize)+1:end,iAnimals)]; %first 50ms shouldnt be used due to deadtime after onset pulse
    end
    try fprintf('done: %s\n',cBhv.Animals{iAnimals}); end
end
