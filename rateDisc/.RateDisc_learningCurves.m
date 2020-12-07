function [Performance,bhv,allDates] = RateDisc_learningCurves(Animal,cPath,lSessions,binSize,highDetection,showPlot,minDelay)
% Analze behavioral data from delayed localization paradigm to get basic readout of animals discrimination performance.
% Reads all available files in 'cPath', ignoring data files that contain less then 100 trials.
%
% Optional inputs:
% lSessions: Only use last 'lSessions' sessions for behavioral analysis.
% binSize: Number of bins used to compute discrimination. Default is 10.
% highDetection: Only use sessions at which detection was at 90% or higher.
% showPlot: Plot behavioral results. Default is true.
% minDelay: Only use sessions at which at least one decisionGap > 'minDelay' was presented. The default is all trials.

%% check optional input
if ~exist('lSessions','var') || isempty(lSessions)
    lSessions = inf;
end

if ~exist('highDetection','var') || isempty(highDetection)
    highDetection = 0;
end

if ~exist('showPlot','var') || isempty(showPlot)
    showPlot = true;
end

if ~exist('minDelay','var') || isempty(minDelay)
    minDelay = 0;
end

if ~exist('binSize','var') || isempty(binSize)
    binSize = 10;
end

[Performance,bhv] = DelayedLoc_learningCurves(Animal,cPath,lSessions,binSize,highDetection,showPlot,minDelay);
Performance.Detection = fliplr(Performance.Detection);
Performance.Discrimination = fliplr(Performance.Discrimination);
Performance.Date = fliplr(Performance.Date);
Performance.DetTrials = fliplr(Performance.DetTrials);
Performance.DiscTrials = fliplr(Performance.DiscTrials);

%% find detection learning curve
prcRange = [5 50 95]; %percentiles to determine learning range. default is 5, 50 and 95th.
detectOn = find(~isnan(Performance.Detection(1,:)),1); %start of detection phase.
discOn = find((~isnan(Performance.Discrimination(1,:))),1) - 2; %start of discrimination phase. Detection is learned at this point.
cDates = Performance.Date(1, detectOn : discOn); %dates for detection sessions
detectCurve = [0.5 Performance.Detection(1, detectOn : discOn)]; %performance in detection sessions
nTrials = [50000 Performance.DetTrials(1, detectOn : discOn)]; %trials in detection sessions

normDates = removeGaps([datenum(cDates(1)); datenum(cDates)] - datenum(cDates(1)) + 1)';
[params, ~, stats] = Behavior_fitPalamedes(normDates ./ max(normDates), detectCurve.*nTrials, nTrials, false, true);

fitRange = 2/max(normDates) : 0.001 : 1; %fit range between first and last detection session
cFit = stats.estimates(3) + (1 - stats.estimates(3) - stats.estimates(4)).*.5.*erfc(-stats.estimates(2).*(fitRange-stats.estimates(1))./sqrt(2));

fitRange = (maxnorm(fitRange) .* (max(normDates)-1))+1;
learnRange = params.guessRate + (1 - params.lapseRate - params.guessRate) * prcRange/100; % get 5, 50 and 95th percentile
[~,a] = max(round(cFit,4)); %find max to restrict percentiles to the dynamic part of the curve
learnLocs = interp1(cFit(1:a), fitRange(1:a), learnRange); %find location of percentiles on the curve

% keep data for plot
Performance.audioLearn.vals = detectCurve(2 : end);
Performance.audioLearn.cFit = cFit;
Performance.audioLearn.fitRange = fitRange;
Performance.audioLearn.learnLocs = learnLocs;
Performance.audioLearn.dates = normDates(2:end);

if showPlot
figure; %make figure for estimates
plot(normDates(2:end),detectCurve(2 : end), 'o', 'linewidth', 2, 'MarkerSize', 8); hold on; ax = gca;
line(ax, [learnLocs(1) learnLocs(1)], [0 1], 'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
line(ax, [learnLocs(2) learnLocs(2)], [0 1], 'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
line(ax, [learnLocs(3) learnLocs(3)], [0 1], 'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
plot(fitRange, cFit, 'linewidth', 4); 
ax.XLim = [0 max(normDates)+1]; axis square; 
ax.YLim = [0.4 1];
xlabel('Sessions'); ylabel('Hit rate');
allDates = {[Animal '; Audio detection'] cDates{round(learnLocs)}};
allDates{2} = [num2str(prcRange(1)) 'th prctile: ' allDates{2}];
allDates{3} = [num2str(prcRange(2)) 'th prctile: ' allDates{3}];
allDates{4} = [num2str(prcRange(3)) 'th prctile: ' allDates{4}];
title(allDates);
end

% add information on expert detection range
allDates{5} = ['Last detection: ' cDates{end}];

%% check for discriminatin learning curve
tactileOn = find((~isnan(Performance.Detection(2,:))),1); %start of somatosensory detection phase.
if isempty(tactileOn)
    tactileOn = size(Performance.Detection,2); %use all discrimination data if no tactile data is found.
end

cDates = Performance.Date(1, discOn + 1 : tactileOn); %dates for audio discrimination sessions
discCurve = Performance.Discrimination(1,  discOn + 1 : tactileOn); %performance in discrimination sessions
nTrials = Performance.DiscTrials(1,  discOn + 1 : tactileOn); %trials in discrimination sessions

normDates = removeGaps(datenum(cDates) - datenum(cDates(1)) + 1)';   
[params, ~, stats] = Behavior_fitPalamedes(normDates ./ max(normDates), discCurve.*nTrials, nTrials, false, true);

fitRange = 1/length(discCurve) : 0.001 : 1; %fit range between first and last detection session
cFit = stats.estimates(3) + (1 - stats.estimates(3) - stats.estimates(4)).*.5.*erfc(-stats.estimates(2).*(fitRange-stats.estimates(1))./sqrt(2));
fitRange = (maxnorm(fitRange) .* (max(normDates)-1))+1;

% keep data for plot
Performance.audioDisc.vals = discCurve;
Performance.audioDisc.cFit = cFit;
Performance.audioDisc.fitRange = fitRange;
Performance.audioDisc.dates = normDates;

if showPlot
figure; %make figure for estimates
plot(normDates,discCurve, 'o', 'linewidth', 2, 'MarkerSize', 8); hold on; ax = gca;
plot(fitRange, cFit, 'linewidth', 4); 
ax.XLim(1) = 1; axis square;
ax.YLim = [0.4 1];
xlabel('Sessions'); ylabel('Hit rate');
allDates{6} = [Animal '; Audio discrimination'];
allDates{7} = ['First audio discrimination: ' cDates{1}];
allDates{8} = ['Last audio discrimination: ' cDates{end}];
title({[Animal ' - Audio discrimination'] allDates{7:8}});
end

%% find tactile detection learning curve
try
tacDetectOn = find(diff(~isnan(Performance.Discrimination(2,:))) == -1, 1); %end of first somatosensory discrimination phase.
allDates{9} = ['Last novice tactile discrimination: ' Performance.Date{tacDetectOn}];

tacDetectOff = find(diff((~isnan(Performance.Discrimination(2,:)))) == 1); %end of somatosensory detection phase.
tacDetectOff = tacDetectOff(tacDetectOff > tacDetectOn + 10);
tacDetectOff = tacDetectOff(1);

cDates = Performance.Date(1, tacDetectOn + 1 : tacDetectOff - 1); %dates fpr detection sessions
detectCurve = [0.5 Performance.Detection(2, tacDetectOn + 1 : tacDetectOff - 1)]; %performance in detection sessions
nTrials = [50000 Performance.DetTrials(2, tacDetectOn + 1 : tacDetectOff - 1)]; %trials in detection sessions

normDates = removeGaps([datenum(cDates(1)); datenum(cDates)] - datenum(cDates(1)) + 1)';
[params, ~, stats] = Behavior_fitPalamedes(normDates ./ max(normDates), detectCurve.*nTrials, nTrials, false, true);

fitRange = 2/length(detectCurve) : 0.001 : 1; %fit range between first and last detection session
cFit = stats.estimates(3) + (1 - stats.estimates(3) - stats.estimates(4)).*.5.*erfc(-stats.estimates(2).*(fitRange-stats.estimates(1))./sqrt(2));

fitRange = (maxnorm(fitRange) .* (max(normDates)-1))+1;
learnRange = params.guessRate + (1 - params.lapseRate - params.guessRate) * prcRange/100; % get 5, 50 and 95th percentile
[~,a] = max(round(cFit,4)); %find max to restrict percentiles to the dynamic part of the curve
learnLocs = interp1(cFit(1:a), fitRange(1:a), learnRange); %find location of percentiles on the curve
if isnan(learnLocs(1)); learnLocs(1) = 1; end
if isnan(learnLocs(3)); learnLocs(3) = fitRange(end); end

% keep data for plot
Performance.tacLearn.vals = detectCurve(2 : end);
Performance.tacLearn.cFit = cFit;
Performance.tacLearn.fitRange = fitRange;
Performance.tacLearn.learnLocs = learnLocs;
Performance.tacLearn.dates = normDates(2:end);

if showPlot
figure; %make figure for estimates
plot(normDates(2:end),detectCurve(2 : end), 'o', 'linewidth', 2, 'MarkerSize', 8); hold on; ax = gca;
plot(fitRange, cFit, 'linewidth', 4); 
line(ax, [learnLocs(1) learnLocs(1)], [0 1], 'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
line(ax, [learnLocs(2) learnLocs(2)], [0 1], 'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
line(ax, [learnLocs(3) learnLocs(3)], [0 1], 'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
ax.XLim = [0 max(normDates)+1]; axis square; 
ax.YLim = [0.4 1];
xlabel('Sessions'); ylabel('Hit rate');
end

allDates{10} = [Animal '; Tactile detection'];
allDates{11} = [num2str(prcRange(1)) 'th prctile: ' cDates{round(learnLocs(1))}];
allDates{12} = [num2str(prcRange(2)) 'th prctile: ' cDates{round(learnLocs(2))}];
allDates{13} = [num2str(prcRange(3)) 'th prctile: ' cDates{round(learnLocs(3))}];
title({[Animal ' - Tactile detection'] allDates{11:13}});

%% check for discriminatin learning curve
tactileDiscOn = find(diff((~isnan(Performance.Discrimination(2,:)))) == 1) + 1; %end of somatosensory detection phase.
tactileDiscOn = min(tactileDiscOn(tactileDiscOn >= tacDetectOff));

tactileDiscOff = find(diff(~isnan(Performance.Discrimination(2,:))) == -1); %end of somatosensory discrimination phase.
tactileDiscOff = min(tactileDiscOff(tactileDiscOff > tactileDiscOn + 5));

cDates = Performance.Date(1,tactileDiscOn : tactileDiscOff); %dates for detection sessions
discCurve = Performance.Discrimination(2, tactileDiscOn : tactileDiscOff); %performance in detection sessions
nTrials = Performance.DiscTrials(2, tactileDiscOn : tactileDiscOff); %trials in detection sessions
[params, ~, stats] = Behavior_fitPalamedes((1 : length(discCurve)) ./ length(discCurve) , discCurve.*nTrials, nTrials, false, true);

fitRange = 1/length(discCurve) : 0.001 : 1; %fit range between first and last detection session
cFit = stats.estimates(3) + (1 - stats.estimates(3) - stats.estimates(4)).*.5.*erfc(-stats.estimates(2).*(fitRange-stats.estimates(1))./sqrt(2));
fitRange = (fitRange - fitRange(1)) ./ (fitRange(end) - fitRange(1)) .* (length(discCurve) - 1) + 1; %rescale fitrange back to sessionNrs

% keep data for plot
Performance.tacDisc.vals = discCurve;
Performance.tacDisc.cFit = cFit;
Performance.tacDisc.fitRange = fitRange;
Performance.tacDisc.dates = datenum';

if showPlot
figure; %make figure for estimates
plot(discCurve, 'o', 'linewidth', 2, 'MarkerSize', 8); hold on; ax = gca;
plot(fitRange, cFit, 'linewidth', 4); 
ax.XLim(1) = 1; axis square;
ax.YLim = [0.4 1];
xlabel('Sessions'); ylabel('Hit rate');
allDates{14} = [Animal '; Tactile discrimination'];
allDates{15} = ['First tactile discrimination: ' Performance.Date{tactileDiscOn}];
allDates{16} = ['Last tactile discrimination: ' Performance.Date{tactileDiscOff}];
title(allDates(14:16));
end

catch
   allDates(9:16) = cell(1,8); % 
end