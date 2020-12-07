function [V, VError] = makeVStruct(folder, M, screenForStability)
% [V, VError] = makeVStruct(folder, M [, screenForStability])
%
% Build a V struct from the alldata files in 'folder', using the units
% specified by M (an M struct). Contains single-trial info in a form to
% make it easy to compute the Fano factor or VarCE. More or less parallel
% to the M struct, but V(u).cond contains a .spikes field, which is an
% nTrials x nTimes logical array of whether that neuron spiked or not at
% that time on that trial. Stimulus-aligned only.
%
% screenForStability determines whether to only use the trials that are
% assessed as 'stable' via markAlldataInstability (as is usually done with
% the M structs). Default 1.

p.preStim = -300;

p.maxMoveDurFast = 1000;  % ms, for use with no-barrier rats
p.maxMoveDurSlow = 2000;  % ms, for use with barrier rats

p.preMove = -200;
p.postMove = 800;

% Should match alldataMeans
p.instabilitySmooth = 10;
p.instabilityThresh = 0.5;
p.instabilityMaxCut = 0.25;
p.instabilityMaxTrend = 0.5;



if ~exist('screenForStability', 'var')
  p.screenForStability = 1;
else
  p.screenForStability = screenForStability;
end

pError = p;
pError.screenForStability = 0;


files = {M.filename};
uFiles = unique(files);


for f = 1:length(uFiles)
  
  %% Load alldata file, narrow down M to units in this alldata file
  
  fprintf('Loading %s... ', uFiles{f});
  
  try
    alldata = loadOneStruct(fullfile(folder, uFiles{f}));
  catch
    error('makeVStruct:alldataLoadFailed', 'Could not load file: %s', fullfile(folder, uFiles{f}));
  end
  
  theseUnits = find(strcmp(files, uFiles{f}));
  thisM = M(theseUnits);
  
  fprintf('adding %d units\n', length(thisM));
  
  
  
  %% Get successful trials
  
  successes = [alldata.hasNeuralData] == 1 & [alldata.earlyWithdrawal] == 0 & [alldata.success] == 1;
  
  V(theseUnits) = makeOneFileV(alldata(successes), thisM, p);
  
  errors = [alldata.hasNeuralData] == 1 & [alldata.earlyWithdrawal] == 0 & [alldata.success] == 0;
  VError(theseUnits) = makeOneFileV(alldata(errors), thisM, pError);
%   alldata = alldata(successes);
end


function V = makeOneFileV(alldata, thisM, p)

%% Eliminate trials with overly long movement times


%   % Add movement duration to data in the process, to make life easy
%   for tr = 1:length(alldata)
%     alldata(tr).movementDur = 1000 * (alldata(tr).parsedEvents.states.reward(1) - ...
%       alldata(tr).parsedEvents.states.tone_playing2(2));
%   end

% Figure out which movement duration cutoff to use
med = nanmedian([alldata.movementDuration]);
if med < 1000
  maxMoveDur = p.maxMoveDurFast;
else
  maxMoveDur = p.maxMoveDurSlow;
end
alldata = alldata([alldata.movementDuration] <= maxMoveDur);


%% Assess stability

if p.screenForStability
%   if ~isfield(alldata(1).units(1), 'stable')
    alldata = markAlldataInstability(alldata, p.instabilitySmooth, p.instabilityThresh, p.instabilityMaxCut, p.instabilityMaxTrend);
%   end
end


%% Construct V

delayLen = min([alldata.imposedWaitDuration]);

% Datasets sometimes use ms and sometimes seconds
if delayLen > 0 && delayLen <= 2
  delayLen = delayLen * 1000;
end

for u = 1:length(thisM)
  V(u).tetrode = thisM(u).tetrode;
  V(u).cluster = thisM(u).cluster;
  V(u).rating = thisM(u).rating;
  
  V(u).cond = getSpikes(alldata, p.preStim, delayLen, p.preMove, p.postMove, thisM(u).tetrode, thisM(u).cluster, p.screenForStability);
  
  V(u).times.stimTimes = p.preStim:delayLen;
  V(u).times.moveTimes = p.preMove:p.postMove;
  
  V(u).stimLen = delayLen;
  V(u).SNR = thisM(u).SNR;
  V(u).filename = thisM(u).filename;
end





function cond = getSpikes(alldata, preStim, delayLen, preMove, postMove, tetrode, cluster, screenForStability)

modes = [-1 0 1];
modeLabels = {'Auditory', 'Multisensory', 'Visual'};
delayLen = floor(delayLen);


% Find the unit
% tetrodes = alldata(1).clusterInfo(:, 1)';
% clusters = alldata(1).clusterInfo(:, 2)';
tetrodes = [alldata(1).units.tetrodeNumber];
clusters = [alldata(1).units.clusterNumber];

u = find(tetrodes == tetrode & clusters == cluster, 1);

% Take only stable trials
if screenForStability
  stable = arrayfun(@(d) d.units(u).stable, alldata);
  alldata = alldata(stable == 1);
end


% Separate left/right and modality
for lr = 1:2
  for mode = 1:length(modes)
    % Modality data
    cond(lr, mode).modality = modes(mode);
    cond(lr, mode).modeLabel = modeLabels{mode};
    
    
    % Get correct trials
    if lr == 1
      cond(lr, mode).leftRight = 'L';
      data = alldata([alldata.visualOrAuditory] == modes(mode) & strcmp({alldata.correctSideName}, 'L'));
    else
      cond(lr, mode).leftRight = 'R';
      data = alldata([alldata.visualOrAuditory] == modes(mode) & strcmp({alldata.correctSideName}, 'R'));
    end
    
    
    % Save rate for each trial
    cond(lr, mode).nVisualEvents = [data.showVisual] .* [data.nVisualEvents];
    cond(lr, mode).nAuditoryEvents = [data.showAudio] .* [data.nAuditoryEvents];
    
    
    % Pre-allocate
    cond(lr, mode).stimAligned.spikes = false(length(data), delayLen - preStim + 1);
    cond(lr, mode).moveAligned.spikes = false(length(data), postMove - preMove + 1);
    
    % Insert spikes
    for tr = 1:length(data)
      S = data(tr).preStimDelay;
      M = data(tr).timeInCenter;
      spikes = round(data(tr).units(u).spikes);
      
      stimSpikes = spikes(spikes >= preStim + S & spikes <= delayLen + S) - preStim - floor(S) + 1;
      moveSpikes = spikes(spikes >= preMove + M & spikes <= postMove + M) - preMove - floor(M) + 1;
      
      cond(lr, mode).stimAligned.spikes(tr, stimSpikes) = true;
      cond(lr, mode).moveAligned.spikes(tr, moveSpikes) = true;
    end
    
  end
end
