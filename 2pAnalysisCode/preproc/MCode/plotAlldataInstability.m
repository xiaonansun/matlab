function plotAlldataInstability(alldata, unit, maxCut)
% plotAlldataInstability(alldata, unit)
% plotAlldataInstability(alldata, unit, maxCut)
%
% Show plots to evaluate unit stability. Takes an alldata struct that has
% been run though markAlldataInstability, and a unit number to display.
% Optionally, takes the 'maxCut' argument used in markAlldataInstability.
% maxCut defaults to 0.25.
%
% Shows three plots. First, it shows a stabilityPlot() for the original
% data for this unit, with the mean for the center portion overlaid.
% Second, it shows a stabilityPlot() for the clipped data with the
% trendSize trend strength, or a blank plot if the whole unit was chucked.
% Last, it shows the stability index with time, with points exceeding
% threshold (black dashed line) marked with blue dots, and if the ends were
% clipped, the clip points are shown with red dots.

%% Optional argument

if ~exist('maxCut', 'var')
  maxCut = 0.25;
end


%% Stability plot for all trials

alldata = alldata([alldata.hasNeuralData] == 1);

% Initial plot of data mean
[FR, trialIDs] = stabilityPlot(alldata, unit);

% Blue bar showing mean
quarters = floor(length(FR) * [maxCut 1-maxCut]);
midMean = mean(FR(quarters(1):quarters(2)));
plot(trialIDs(floor(length(FR) * [maxCut 1-maxCut])), midMean * [1 1], 'b');

title('Original');


%% Grab the key values

stable = arrayfun(@(x) x.units(unit).stable, alldata);

index = arrayfun(@(x) x.units(unit).stabilityVals.instabilityIndex, alldata);

trendVals = arrayfun(@(x) x.units(unit).stabilityVals.trendSize, alldata);
trendVals = trendVals(~isnan(trendVals));
if ~isempty(trendVals)
  trend = trendVals(1);
else
  trend = NaN;
end

thresh = alldata(1).stabilityParams.instabilityThresh;


%% 'After treatment' stability plot

% Plot trimmed data or generate blank figure to show unit dropped
anyStable = any(stable);
if anyStable
  stabilityPlot(alldata(stable), unit);
  plot(trialIDs(floor(length(FR) * [maxCut 1-maxCut])), midMean * [1 1], 'b');
  title('After cleanup');
  
else
  figure;
  axis;
  title('Unit thrown away');
end
textTitle(sprintf('Trend = %0.2f', trend));


%% Plot instability index

nonNans = find(~isnan(index));
indNoNans = index(nonNans);
overThresh = indNoNans > thresh;

figure;
hold on;
% Plot threshold
plot([1 nonNans(end)], thresh * [1 1], 'k--');
% Plot index
plot(nonNans, indNoNans);
% Plot points above threshold with dots
plot(nonNans(overThresh), indNoNans(overThresh), 'b.');

% Plot clip points in red dots if the unit wasn't chucked
if anyStable
  % Find the first stable point, then the last non-NaN unstable point
  % preceeding it
  firstStable = find(stable, 1);
  prevUnstable = find(~isnan(index(1:firstStable-1)), 1, 'last');
  
  % Same for last stable point
  lastStable = find(stable, 1, 'last');
  follUnstable = lastStable + find(~isnan(index(lastStable+1:end)), 1);
  
  % Plot red dots if clipped
  if  ~isempty(prevUnstable)
    plot(prevUnstable, index(prevUnstable), 'r.');
  end
  if ~isempty(follUnstable)
    plot(follUnstable, index(follUnstable), 'r.');
  end
end

% Set y-axis to go down to zero, so that the scale is reasonable
yLims = get(gca, 'YLim');
set(gca, 'YLim', [0 yLims(2)]);

title('Instability index');

