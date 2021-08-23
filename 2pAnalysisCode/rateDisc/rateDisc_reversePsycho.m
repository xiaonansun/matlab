[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
[cTrainDates, cTrainLabels] = rateDisc_convertTrainDates(trainDates);
trainDur = 'AudioDisc'; %period from behavioral data should be used
cPath = 'Y:\data\Behavior_Simon\';
Animals = {'mSM63' 'mSM64' 'mSM65' 'mSM66'};
binSize = 0.025; %10 ms steps
distBins = 10;

%% collect data
bhv = []; Cnt = 0;
for iAnimals = 1 : length(Animals)
    
    clear recs
    recFiles = dir([cPath Animals{iAnimals} filesep 'SpatialDisc' filesep 'Session Data' filesep]); %behavioral data from current animal
    for iFiles = 1 : length(recFiles)
        a = textscan(recFiles(iFiles).name, '%s%s%s%s%s', 'Delimiter','_');
        if ~isempty(a{3})
            recs(iFiles).name = [a{3}{1} '_' a{4}{1}]; %isolate date of current recording
        end
    end
    
    [recs, useIdx] = rateDisc_filterRecs(trainDates{iAnimals}, recs, trainDur); %selected recordings from designated training range
    recFiles = recFiles(useIdx); % behavioral files that should be used

    for iFiles = 1 : length(recFiles)
        load([cPath Animals{iAnimals} filesep 'SpatialDisc' filesep 'Session Data' filesep recFiles(iFiles).name])
        useData = length(SessionData.Rewarded) > 100 && length(unique(SessionData.DistStim)) > 1 && sum(SessionData.DistStim > 0) > 50; % if file contains at least 100 trials and different distractor rates;
        if useData
            
            Cnt = Cnt + 1;
            SessionData.SessionNr = repmat(Cnt,1,SessionData.nTrials); %tag all trials in current dataset with session nr
            SessionData.AnimalNr = repmat(iAnimals,1,SessionData.nTrials); %tag all trials in current dataset with session nr
            if mean(SessionData.Rewarded(1:100)) < 0.6 %don't use intiial trials if performance is too low
                SessionData = selectBehaviorTrials(SessionData, 101 : size(SessionData.Rewarded,2));
            end
            bhv = appendBehavior(bhv,SessionData); %append into larger array
        end
    end
end

%% run logistic regresion
clear lWeights rWeights
for iAnimals = 1 : length(Animals) + 1
    
    % stimDur = round(max(bhv.stimDur)/binSize) * binSize; %make sure stimDur is can be divided by binSize
    stimDur = 1;
    leftRate = zeros(length(bhv.stimEvents),stimDur/binSize); %allocate leftRate array
    rightRate = zeros(length(bhv.stimEvents),stimDur/binSize); %allocate rightRate array
    
    for iTrials = 1:length(bhv.stimEvents)
        leftRate(iTrials,:) = histcounts(bhv.stimEvents{iTrials}{1},0:binSize:stimDur);  %count events in each bin for left events
        rightRate(iTrials,:) = histcounts(bhv.stimEvents{iTrials}{2},0:binSize:stimDur);  %count events in each bin for right events
    end
    
    excessRate = leftRate-rightRate; %compute excess rate
    distFractions = bhv.DistStim ./ bhv.TargStim;
    lInd = bhv.Assisted & ((bhv.CorrectSide == 1 & bhv.Rewarded) |(bhv.CorrectSide == 2 & bhv.Punished)); %index for non-detection, left-choice trials
    rInd = bhv.Assisted & ((bhv.CorrectSide == 2 & bhv.Rewarded) |(bhv.CorrectSide == 1 & bhv.Punished)); %index for non-detection, right-choice trials
    
    usedRate = 0;
    if iAnimals > length(Animals)
        trialSelect = distFractions >= usedRate;
    else
        trialSelect = distFractions >= usedRate & bhv.AnimalNr == iAnimals;
    end
    lWeights(:, iAnimals) = mnrfit(excessRate(trialSelect,:),lInd(trialSelect)'+1);
    rWeights(:, iAnimals) = mnrfit(excessRate(trialSelect,:),rInd(trialSelect)'+1);
    lWeights(:, iAnimals) = [zeros(round(0.05/binSize),1);lWeights(round(0.05/binSize)+1:end,iAnimals)]; %first 50ms shouldnt be used due to deadtime after onset pulse
    rWeights(:, iAnimals) = [zeros(round(0.05/binSize),1);rWeights(round(0.05/binSize)+1:end,iAnimals)]; %first 50ms shouldnt be used due to deadtime after onset pulse
    
end

%%
% plot(binSize:binSize:length(lWeights)*binSize,lWeights,'g','linewidth',2);hold on
figure;
stdshade(lWeights',0:binSize:stimDur,'g',0.5,3); hold on
stdshade(rWeights',0:binSize:stimDur,'r',0.5,3);
title(['Logistic regression weights']);
axis square; hline(0); 
ylabel('beta weights');xlabel('time(s)');
   
%% compute discrimination performance
targMod = 2; %target modality (1 = vision, 2 = audio, 4 = somatosensory)
targAnimal = 4;
ind = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.StimType == targMod & bhv.AnimalNr == targAnimal; %only use active trials
rInd = (bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded); %right-choice trials
rInd = rInd(ind);

tCnt = 0;
eventCnt = zeros(2,sum(ind));
for iTrials = find(ind)
    tCnt = tCnt +1;
    [left, right] = Behavior_getStimEvent(bhv.StimType(iTrials), bhv.stimEvents{iTrials});
    eventCnt(1,tCnt) = length(right);
    eventCnt(2,tCnt) = length(left) + length(right);
end

[nTrials, distRatio, trialInd]  = histcounts(eventCnt(1,:) ./ eventCnt(2,:), distBins); %get ratio between rightward and absolute number of events
distRatio = distRatio + diff(distRatio(1:2))/2; distRatio(end) = [];
for iBins = 1:length(nTrials)
    rightChoice(iBins) = sum(rInd(trialInd == iBins)); %get number of rightward trials for each difficulty
end
[params, h1, ~, cFit] = Behavior_fitPalamedes(distRatio, rightChoice, nTrials, true, true); %complete discrimination parameters
title(h1.data.Parent,[Animals{targAnimal} ' - All perf. trials'])
ylim(h1.data.Parent, [0 1]);