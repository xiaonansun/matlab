function hf = decisionPSTH(M, showErrBars, titl, plotPerModality, spawnNewFigure, alphaVal)
% hf = decisionPSTH(M)
% hf = decisionPSTH(M, showErrBars)
% hf = decisionPSTH(M, showErrBars, titl)
% hf = decisionPSTH(M, showErrBars, titl, plotPerModality)
% hf = decisionPSTH(M, showErrBars, titl, plotPerModality, spawnNewFigure)
% hf = decisionPSTH(M, showErrBars, titl, plotPerModality, spawnNewFigure, alphaVal)
%
% Plot a nice PSTH of a unit, specified by an M struct. hf is the handle to
% the figure plotted in.
%
% Optional arguments:
%
% showErrBars, if 1, displays errorbars; if -1, displays error fills; if 0,
%   no errors shown. Default -1.
% titl, if a non-empty string, puts that title on the plot. If missing or
%   empty, the title shows the tetrode, cluster, rating, and SNR. To have
%   no title, make the string, for example, ' '.
% plotPerModality, if 1, plots one PSTH per modality. Default 0.
% spawnNewFigure, if 1, spawns a new figure before plotting. Default 1.
% alphaVal determines the alpha transparency for error fills. Default 0.15.
%   Use alphaVal = 1 to enable saving as an Illustrator file.

%% Parameters
modalityColors = {[0 1 0], [0 0 1], [0 0 0]};
lrStyles = {'-', '--'};

if ~exist('showErrBars', 'var')
  showErrBars = -1;
end

if ~exist('titl', 'var') || isempty(titl)
  if ~isfield(M, 'filename')
    titl = sprintf('Tetrode %d, cluster %d, rating %0.1f, SNR %1.1f', M.tetrode, M.cluster, M.rating, M.SNR);
  else
    fixedFile = regexprep(M.filename, '_', '\\_');
    titl = sprintf('Tetrode %d, cluster %d, rating %0.1f, SNR %1.1f\n%s', M.tetrode, M.cluster, M.rating, M.SNR, fixedFile);
  end
end

if ~exist('plotPerModality', 'var')
  plotPerModality = 0;
end

if ~exist('spawnNewFigure', 'var')
  spawnNewFigure = 1;
end

if ~exist('alphaVal', 'var')
  alphaVal = 0.15;
end


%% Error checking
if length(M) > 1
  error('Provide only one neuron');
end


%% If no move data, mark times so that move-locked axis won't display
if isempty(M.cond(1).moveAligned.meanFR)
  M.times.moveTimes(end) = -Inf;
end


%% Loop through modalities

nModes = size(M.cond, 2);

if spawnNewFigure
  hf = blankFigure;
else
  hf = gca;
end

hold on;

theMax = 0;
for mode = 1:nModes
  % If needed, start a new figure
  if plotPerModality ~= 0 && mode > 1
    theMax = 0;
    
    if spawnNewFigure
      hf(end+1) = blankFigure;
    else
      hf(end+1) = gca;
    end
    
    hold on;
  end
  
  %% Loop through strengths
  
  nStrengths = size(M.cond(1, mode, :), 3);
  
  for s = 1:nStrengths
    color = modalityColors{mode} * (1 + nStrengths - s) / nStrengths + [1 1 1] * (s - 1) / nStrengths;
    
    for lr = 1:2
      thisMax = plotLines(M.cond(lr, mode, s), M.times, showErrBars, ...
        color, lrStyles{lr}, 2, alphaVal);
      theMax = max(theMax, thisMax);
    end
  end
  
  if plotPerModality
    addAxes(M.times, theMax);
    textTitle([titl, ' ', M.cond(1, mode).modeLabel]);
%     addRTDist(M.cond(1:2, mode), M.times, newMax);
  end
end

if ~plotPerModality
  newMax = addAxes(M.times, theMax);
  if spawnNewFigure
    textTitle(titl);
    xLims = get(gca, 'XLim');
    y = newMax * 0.95;
    for m = 1:nModes
      text(xLims(2), y, M.cond(1, m).modeLabel, 'color', modalityColors{m}, 'HorizontalAlignment', 'right');
      y = y - newMax * 0.05;
    end
  end
%   addRTDist(M.cond, M.times, newMax);
end


function theMax = plotLines(vals, times, showErrBars, color, style, mainThick, alphaVal)

spacing = 200;

% Stim-aligned
plot(times.stimTimes, vals.stimAligned.meanFR, 'color', color, 'lineStyle', style, 'lineWidth', mainThick);

% Deal with spacing
lastStimTime = times.stimTimes(end);
firstMoveTime = times.moveTimes(1);
spacingAdd = lastStimTime - firstMoveTime + spacing;

% Move-aligned
if ~isempty(vals.moveAligned.meanFR)
  plot(times.moveTimes + spacingAdd, vals.moveAligned.meanFR, ...
    'color', color, 'lineStyle', style, 'lineWidth', mainThick);
end


if showErrBars == 1
  plot(times.stimTimes, vals.stimAligned.meanFR + vals.stimAligned.stdErr, ...
    'color', color, 'lineStyle', style, 'lineWidth', 0.5);
  plot(times.stimTimes, vals.stimAligned.meanFR - vals.stimAligned.stdErr, ...
    'color', color, 'lineStyle', style, 'lineWidth', 0.5);
  
  if ~isempty(vals.moveAligned.meanFR)
    plot(times.moveTimes + spacingAdd, vals.moveAligned.meanFR + vals.moveAligned.stdErr, ...
      'color', color, 'lineStyle', style, 'lineWidth', 0.5);
    plot(times.moveTimes + spacingAdd, vals.moveAligned.meanFR - vals.moveAligned.stdErr, ...
      'color', color, 'lineStyle', style, 'lineWidth', 0.5);
  end
  
elseif showErrBars == -1
  % Error fills
  
  sTimes = [times.stimTimes fliplr(times.stimTimes)];
  eVals = [vals.stimAligned.meanFR + vals.stimAligned.stdErr, ...
    fliplr(vals.stimAligned.meanFR - vals.stimAligned.stdErr)];
  h = fill(sTimes, eVals, color, 'LineStyle', 'none');
  if alphaVal < 1
    alpha(h, alphaVal);
  end
  
  
  if ~isempty(vals.moveAligned.meanFR)
    mTimes = [times.moveTimes fliplr(times.moveTimes)] + spacingAdd;
    eVals = [vals.moveAligned.meanFR + vals.moveAligned.stdErr, ...
      fliplr(vals.moveAligned.meanFR - vals.moveAligned.stdErr)];
    
    h = fill(mTimes, eVals, color, 'LineStyle', 'none');
    alpha(h, alphaVal);
  end
end


if ~isempty(vals.moveAligned.meanFR)
  theMax = max(max(vals.stimAligned.meanFR), max(vals.moveAligned.meanFR));
else
  theMax = max(vals.stimAligned.meanFR);
end



function theMax = addAxes(times, theMax)

spacing = 200;  % Must match spacing in plotLines above

%% Y axis

theMax = ceil(theMax);

if theMax < 5
  tickInt = theMax;
elseif theMax < 10
  tickInt = 5;
elseif theMax < 30
  tickInt = 10;
elseif theMax < 60
  tickInt = 20;
else
  tickInt = 25;
end

theMax = tickInt * ceil(theMax / tickInt);

ticks = 0:tickInt:theMax;

params.axisOrientation = 'v';
params.fontSize = 12;
params.tickLocations = ticks;
params.axisLabel = 'spikes/s';
params.axisOffset = times.stimTimes(1) - 50;

AxisMMC(0, theMax, params);


%% X-axes

clear params

% Stim-locked
params.axisOrientation = 'h';
params.fontSize = 12;
params.axisOffset = -theMax / 50;
params.tickLocations = [times.stimTimes(1) 0 times.stimTimes(end)];
params.tickLength = 0.015 * theMax;
params.tickLabelLocations = params.tickLocations;
params.tickLabels = {num2str(times.stimTimes(1)) 'Stim' num2str(times.stimTimes(end))};

AxisMMC(params.tickLocations(1), params.tickLocations(end), params);

% Move-locked

if ~isinf(times.moveTimes(end))
  % Deal with spacing
  lastStimTime = times.stimTimes(end);
  firstMoveTime = times.moveTimes(1);
  spacingAdd = lastStimTime - firstMoveTime + spacing;
  
  params.tickLocations = spacingAdd + [times.moveTimes(1) 0 times.moveTimes(end)];
  params.tickLabelLocations = params.tickLocations;
  params.tickLabels = {num2str(times.moveTimes(1)) 'Move' num2str(times.moveTimes(end))};
  
  AxisMMC(params.tickLocations(1), params.tickLocations(end), params);
end
