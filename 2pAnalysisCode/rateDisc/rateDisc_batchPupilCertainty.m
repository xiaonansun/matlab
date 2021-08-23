[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
animals = dataOverview(1:8,1);
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
baseDur = 2;          % Duration of baseline before lever grab in seconds
sRate = 30;

bhv = [];
allPupil = [];
allTimes = [];
for iAnimals = 1 : length(animals)
    
    cAnimal = animals{iAnimals};
    recs = dir([cPath filesep cAnimal filesep 'SpatialDisc' filesep]);
    recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs, 'audioDisc'); %this sorts recordings by date
    fprintf('Basepath: %s; Found %d recordings\n', cPath, length(recs));
    
    for iRecs = 1 : length(recs)
        try
            disp([cAnimal ' - ' recs(iRecs).name]);
            fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %session data path
            
            [cBhv, cPupil, cTimes] = rateDisc_pupilCertainty(fPath, baseDur);
            cBhv.AnimalID = ones(1, length(cBhv.Rewarded)) * iAnimals;
            cBhv.SessionNr = ones(1, length(cBhv.Rewarded)) * iRecs;
            bhv = appendBehavior(bhv,cBhv); %append into larger array
            
            allPupil = [allPupil cPupil];
            allTimes = [allTimes cTimes];
            
        end
    end
end

%% get previous choice bias for discrimination
distBins = 6;
% dInd = ismember(bhv.AnimalID,1:8);
dInd = ismember(bhv.AnimalID, [1 2 3 4]);

figure
subplot(1,2,1)
%Plot fit and real data
[distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, dInd, distBins);
plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', 'k'); hold on
errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', 'k', 'MarkerFaceColor','w','linewidth',2);

% rInd = dInd & bhv.CorrectSide == 2; %right-choice trials
rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
rInd = [rInd(2:end) false];
[distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, rInd, distBins);
plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', 'g'); hold on
errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', 'g', 'MarkerFaceColor','w','linewidth',2);

% lInd = dInd & bhv.CorrectSide == 1; %right-choice trials
lInd = dInd & ~((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
lInd = [lInd(2:end) false];
[distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, lInd, distBins);
plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', 'r'); hold on
h = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', 'r', 'MarkerFaceColor','w','linewidth',2);

xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
ylim([0 1]); ylabel('Proportion chose right'); hold off; axis square
xlabel('Distractor ratio'); title('Frontal bilateral inactivation');
disp(nTrials);
niceFigure(h.Parent);

% get previous choise bias in general
% dInd = ismember(bhv.AnimalID,1:8);
dInd = ismember(bhv.AnimalID, [1 2 3 4]);

subplot(1,2,2);
rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
lInd = dInd & ~((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
prInd = [rInd(2:end) false];
plInd = [lInd(2:end) false];

% all trials
[rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd), sum(dInd));

% previous trial
[prHigh, prLow, prR] = Behavior_wilsonError(sum(rInd(prInd)), sum(prInd));
[plHigh, plrLow, plR] = Behavior_wilsonError(sum(rInd(plInd)), sum(plInd));

data = [prR rProb plR];
dataHigh = [prHigh rHigh plHigh];
dataLow = [prLow rLow plrLow];

h = errorbar(1:length(data), data, data-dataLow, dataHigh-data, 'o', 'color', 'k', 'MarkerFaceColor','w','linewidth',2);
niceFigure(h.Parent); axis square; xlim([0 length(data)+1]); hline(0.5);
title('Previous choice bias');

%%
dInd = ismember(bhv.AnimalID,1:8);
rInd = dInd & ((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
lInd = dInd & ~((bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded)); %right-choice trials
prInd = [rInd(2:end) false];
plInd = [lInd(2:end) false];

figure; subplot(1,2,1); hold on;
Cnt = 0; clear cBias cHigh cLow
for iBins = 1 : size(allPupil,1) - 5
    Cnt = Cnt + 1;
    cPupil = nanmean(allPupil(iBins : iBins + 5, :), 1);
    
    
    lowInd = cPupil <= prctile(cPupil(dInd), 33) & dInd;
    highInd = cPupil >= prctile(cPupil(dInd), 66) & dInd;
    
    [rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd(lowInd)), sum(lowInd));
    [cHigh(Cnt,1), cLow(Cnt,1), cBias(Cnt,1)] = Behavior_wilsonError((sum(rInd(prInd & lowInd)) + sum(lInd(plInd & lowInd))), sum(lowInd));
    [cHigh(Cnt,2), cLow(Cnt,2), cBias(Cnt,2)] = Behavior_wilsonError((sum(rInd(prInd & highInd)) + sum(lInd(plInd & highInd))), sum(highInd));
    
end
useIdx = ~(isnan(cBias(:,1)) | isnan(cBias(:,2)));
F = (1 : sum(useIdx));
fill([F fliplr(F)],[cHigh(useIdx,1)' fliplr(cLow(useIdx,1)')], 'g', 'linestyle','none', 'FaceAlpha', 0.5);
plot(cBias(useIdx,1), 'linewidth', 4, 'color', 'g');

fill([F fliplr(F)],[cHigh(useIdx,2)' fliplr(cLow(useIdx,2)')], 'b', 'linestyle','none', 'FaceAlpha', 0.5);
plot(cBias(useIdx,2), 'linewidth', 4, 'color', 'b');

%%
baseWin = 5.5;
allWins = [3 1 1];

figure; 
subplot(1,4,1); hold on;
Cnt = 0; clear cBias cHigh cLow
for iBins = 1 : baseWin * sRate
    try
    Cnt = Cnt + 1;
    cPupil = nanmean(allPupil(iBins : iBins + 5, :), 1);
    cPupil = [cPupil(2:end) NaN];
    
    lowInd = cPupil <= prctile(cPupil, 33);
    highInd = cPupil >= prctile(cPupil, 66);
    
    [rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd(lowInd)), sum(lowInd));
    [cHigh(Cnt,1), cLow(Cnt,1), cBias(Cnt,1)] = Behavior_wilsonError((sum(rInd(prInd & lowInd)) + sum(lInd(plInd & lowInd))), sum(lowInd));
    [cHigh(Cnt,2), cLow(Cnt,2), cBias(Cnt,2)] = Behavior_wilsonError((sum(rInd(prInd & highInd)) + sum(lInd(plInd & highInd))), sum(highInd));
    end
end
useIdx = ~(isnan(cBias(:,1)) | isnan(cBias(:,2)));
F = (1 : sum(useIdx));
fill([F fliplr(F)],[cHigh(useIdx,1)' fliplr(cLow(useIdx,1)')], 'g', 'linestyle','none', 'FaceAlpha', 0.5);
plot(cBias(useIdx,1), 'linewidth', 4, 'color', 'g');

fill([F fliplr(F)],[cHigh(useIdx,2)' fliplr(cLow(useIdx,2)')], 'b', 'linestyle','none', 'FaceAlpha', 0.5);
plot(cBias(useIdx,2), 'linewidth', 4, 'color', 'b');
niceFigure(gca); axis square;

% pupil effects for different trial windows
for x = 1 : length(allWins)
    subplot(1,4,x+1); hold on;
    Cnt = 0; clear cBias cHigh cLow
    for iBins = -29 : allWins(x) * sRate
        try
        Cnt = Cnt + 1;
        cPupil = nanmean(allPupil(allTimes(x,:) + iBins : allTimes(x,:) + iBins + 5, :), 1);
        cPupil = [cPupil(2:end) NaN];
        
        lowInd = cPupil <= prctile(cPupil, 33);
        highInd = cPupil >= prctile(cPupil, 66);
        
        [rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd(lowInd)), sum(lowInd));
        [cHigh(Cnt,1), cLow(Cnt,1), cBias(Cnt,1)] = Behavior_wilsonError((sum(rInd(prInd & lowInd)) + sum(lInd(plInd & lowInd))), sum(lowInd));
        [cHigh(Cnt,2), cLow(Cnt,2), cBias(Cnt,2)] = Behavior_wilsonError((sum(rInd(prInd & highInd)) + sum(lInd(plInd & highInd))), sum(highInd));
    end
    end
    useIdx = ~(isnan(cBias(:,1)) | isnan(cBias(:,2)));
    F = (1 : sum(useIdx));
    fill([F fliplr(F)],[cHigh(useIdx,1)' fliplr(cLow(useIdx,1)')], 'g', 'linestyle','none', 'FaceAlpha', 0.5);
    plot(cBias(useIdx,1), 'linewidth', 4, 'color', 'g');
    
    fill([F fliplr(F)],[cHigh(useIdx,2)' fliplr(cLow(useIdx,2)')], 'b', 'linestyle','none', 'FaceAlpha', 0.5);
    plot(cBias(useIdx,2), 'linewidth', 4, 'color', 'b');
    niceFigure(gca); axis square;
end

%%
baseWin = 5.5;
allWins = [2 1 1];

figure; 
subplot(1,4,1); hold on;
Cnt = 0; clear cBias cHigh cLow
for iBins = 1 : baseWin * sRate
    try
    Cnt = Cnt + 1;
    cPupil = nanmean(allPupil(iBins : iBins + 5, :), 1);
    cPupil = [cPupil(2:end) NaN];
    
    lowInd = cPupil <= prctile(cPupil, 33);
    highInd = cPupil >= prctile(cPupil, 66);
    
    [rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd(lowInd)), sum(lowInd));
    [cHigh(Cnt,1), cLow(Cnt,1), cBias(Cnt,1)] = Behavior_wilsonError((sum(rInd(prInd & lowInd)) + sum(lInd(plInd & lowInd))), sum(lowInd));
    [cHigh(Cnt,2), cLow(Cnt,2), cBias(Cnt,2)] = Behavior_wilsonError((sum(rInd(prInd & highInd)) + sum(lInd(plInd & highInd))), sum(highInd));
    end
end
useIdx = ~(isnan(cBias(:,1)) | isnan(cBias(:,2)));
F = (1 : sum(useIdx));
fill([F fliplr(F)],[cHigh(useIdx,1)' fliplr(cLow(useIdx,1)')], 'g', 'linestyle','none', 'FaceAlpha', 0.5);
plot(cBias(useIdx,1), 'linewidth', 4, 'color', 'g');

fill([F fliplr(F)],[cHigh(useIdx,2)' fliplr(cLow(useIdx,2)')], 'b', 'linestyle','none', 'FaceAlpha', 0.5);
plot(cBias(useIdx,2), 'linewidth', 4, 'color', 'b');
niceFigure(gca); axis square;

% pupil effects for different trial windows
a = cumsum(repmat(size(allPupil,1),1,size(allPupil,2))) - size(allPupil,1);
for x = 1 : length(allWins)
    subplot(1,4,x+1); hold on;
    Cnt = 0; clear cBias cHigh cLow
    for iBins = -30 : allWins(x) * sRate
        
        cInd = repmat(allTimes(1,:) + a + iBins,5,1) + (0:4)';
        
        cPupil = nanmean(allPupil(cInd(:)), 1);
        
        
        try
            Cnt = Cnt + 1;
            %             cPupil = nanmean(allPupil(allTimes(x,:) + iBins : allTimes(x,:) + iBins + 5, :), 1);
            
            cAvg = [];
            for xx = 1 :5
                cPupil = allPupil(allTimes(1,:) + a + iBins + xx);
                cAvg = runMean(cAvg, cPupil, xx);
            end
            cPupil = [cAvg(2:end) NaN];
            
            lowInd = cPupil <= prctile(cPupil, 33);
            highInd = cPupil >= prctile(cPupil, 66);
            
            [rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd(lowInd)), sum(lowInd));
            [cHigh(Cnt,1), cLow(Cnt,1), cBias(Cnt,1)] = Behavior_wilsonError((sum(rInd(prInd & lowInd)) + sum(lInd(plInd & lowInd))), sum(lowInd));
            [cHigh(Cnt,2), cLow(Cnt,2), cBias(Cnt,2)] = Behavior_wilsonError((sum(rInd(prInd & highInd)) + sum(lInd(plInd & highInd))), sum(highInd));
        end
    end
    useIdx = ~(isnan(cBias(:,1)) | isnan(cBias(:,2)));
    F = (1 : sum(useIdx));
    fill([F fliplr(F)],[cHigh(useIdx,1)' fliplr(cLow(useIdx,1)')], 'g', 'linestyle','none', 'FaceAlpha', 0.5);
    plot(cBias(useIdx,1), 'linewidth', 4, 'color', 'g');
    
    fill([F fliplr(F)],[cHigh(useIdx,2)' fliplr(cLow(useIdx,2)')], 'b', 'linestyle','none', 'FaceAlpha', 0.5);
    plot(cBias(useIdx,2), 'linewidth', 4, 'color', 'b');
    niceFigure(gca); axis square;
end

%% make error bar plot for different trial episodes
clear cIdx cHigh cLow cBias
cIdx{1} = 60 : 105; %stim
cIdx{2} = 106 : 135; %delay
cIdx{3} = 136 : 160; %response

for iBins = 1 : length(cIdx)

cPupil = nanmean(allPupil(cIdx{iBins}, :), 1);
cPupil = [cPupil(2:end) NaN];

lowInd = cPupil <= prctile(cPupil, 33);
highInd = cPupil >= prctile(cPupil, 66);
    
[rHigh, rLow, rProb] = Behavior_wilsonError(sum(rInd(lowInd)), sum(lowInd));
[cHigh(iBins,1), cLow(iBins,1), cBias(iBins,1)] = Behavior_wilsonError((sum(rInd(prInd & lowInd)) + sum(lInd(plInd & lowInd))), sum(lowInd));
[cHigh(iBins,2), cLow(iBins,2), cBias(iBins,2)] = Behavior_wilsonError((sum(rInd(prInd & highInd)) + sum(lInd(plInd & highInd))), sum(highInd));

end

figure; hold on;
bar(cBias);
errorbar((1:iBins) - 0.15, cBias(:,1), cBias(:,1) - cLow(:,1), cHigh(:,1) - cBias(:,1), 'ko');
errorbar((1:iBins) + 0.15, cBias(:,2), cBias(:,2) - cLow(:,2), cHigh(:,2) - cBias(:,2), 'ko');
ax = gca; ax.XTick = 1:iBins; 
ax.XTickLabel = {'Stimulus' 'Delay' 'Response'};
niceFigure(gca); ylim([0.5 0.6]); axis square
