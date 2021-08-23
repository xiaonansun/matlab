function data = twoP_alignDetectionTask(ops, spks, iscell, redcell, binFilename, showPlots)
% data = twoP_alignDetectionTask(datFilename, binFilename [, showPlots])
% data = twoP_alignDetectionTask(dat, binFilename [, showPlots])
% 

% 2020-04-24 This script is updated from a previous version written by Simon
% Musall and Matt Kaufman for compatiblity with the newest Python version of
% Suite2p. Revisions include: (1) dFOF is no no longer a not part of the
% data struct, (2) S is the spks.npy output from suite2p, (3) DS is no long
% part of the suite2p output, it is instead concatenated from ops.xoff and
% ops.yoff.

% For a 2P session in Simon's detection task: load the analog channels
% (binFilename) and the neural data output from Suite2P (datFilename), and
% convert to a trial-aligned format. If you've already loaded the neural
% data, you may supply it instead of the filename.
% 
% Optional argument showPlots shows some sanity-check plots. Default false.
% 
% Output structure data has fields:
%     analogIdentities: names of the analog channels
%           analogFreq: sampling frequency of the analog data
%         trialNumbers: for matching with behavioral data
%                    A: neuron spatial footprints (pixels x neurons)
%                   im: projection image of the data (maybe mean?)
%       trialStimFrame: the frame of the aligned neural data when the
%                       stimulus occurred
%               neural: aligned neural data. neurons x times x trials
%          neuralTimes: times of the neural data relative to stimulus, in
%                       frames
%        neuralTimesMs: times of the neural data relative to stimulus, in ms
%      trialStimSample: the sample of the aligned analog data when the
%                       stimulus occurred
%               analog: aligned analog data. channels x times x trials
%          analogTimes: times of the analog data relative to stimulus, in
%                       samples (which is also ms)
%                 dFOF: deltaF/F of neural data.
%                   DS: Deviation from reference frame in x and y in each
%                       frame. Indicates 2D motion.
% 
% Parts of this code based on framesPerTrialStopStart3An.

%% commentable variables
% ops = npy.ops;
% spks = npy.spks;
% iscell = npy.iscell;
% redcell = npy.redcell;
% binFilename = npy.bin_MScan_filepath;


%% Parameters

% Sampling

% This value is shown (as frame duration) in the File Properties button of
% File Manager window of Mview, presumably more accurate than 1000/30.9 Hz
msPerFrame =  32.363833;

data.msPerFrame = msPerFrame;

% Time window
preTime = 3;
postTime = 4;

% TTL interpretation
voltageThresh = 1.5;  % V
% Simon used longer intervals for barcodes than Matt Kaufman's defaults
shortInt = 2.5;
longInt = 5.5;

% Analog channel identities

trialCodeCh = 1;
slowGalvoCh = 2;
trialOffCh = 3;
stimOnCh = 4;

analogChannels = 1:5;  % channels to save
data.analogIdentities = {'trialCode', 'slowGalvo', 'trialOff', 'stimOn', 'piezo'};

data.animal = ops.animal;
data.session = ops.session;

%% Optional arguments

if ~exist('showPlots', 'var')
  showPlots = 0;
end


%% Load in analog data

tic
[volt, data.analogFreq] = readMOMAnalog(binFilename);
disp(['MScan analog file (.bin) loaded in ' num2str(toc) ' seconds.']);

trialCodes = volt(trialCodeCh, :);
slowGalvo = volt(slowGalvoCh, :);
trialStop = volt(trialOffCh, :);
stimOn = volt(stimOnCh, :);


%% Find alignment values, decode trial codes

stimOn(stimOn <= voltageThresh) = 0;
stimOn(stimOn > voltageThresh) = 1;
stimOnSamples = find(diff(stimOn) == 1) + 1;

[trialNumbers, codeStarts] = segmentVoltageAndReadBarcodes(trialCodes, shortInt, longInt);

trialStop(trialStop < voltageThresh) = 0;
trialStop(trialStop > voltageThresh) = 1;
trialStopSamples = find(diff(trialStop) == 1) + 1;

%% check if trialNumbers are monotonic. If not if they could belong to multiple sessions.
trials = trialNumbers(trialNumbers > -1);
nrChange = find(diff(trials) ~= 1);

if ~isempty(nrChange)
    if nrChange > 100 || length(nrChange) > 1
        warning('Trial numbers are not monotonic! Increment subsequent trialnumbers to recover full session!')
        for iChange = 1 : length(nrChange)
            trials(nrChange + 1 : end) = trials(nrChange + 1 : end) + trials(nrChange);
        end
        trialNumbers(trialNumbers > -1) = trials;
    else
        if length(trials) == length(stimOnSamples)
            warning('Trial numbers are not monotonic after less then 100 trials. Rejecting those trials!');
            stimOnSamples(1:nrChange) = [];
        else
            error('Trial numbers are not monotonic and dont match stimulus onset times');
        end
    end
end


%% Now, find the starts, stops, and trial codes corresponding to each stimOn

[starts, stops, data.trialNumbers] = alignStartsStopsNumbers(stimOnSamples, codeStarts, trialNumbers, trialStopSamples);
if any(isnan(starts))
    stops(isnan(starts)) = [];
    data.trialNumbers(isnan(starts)) = [];
    starts(isnan(starts)) = [];
end

%% Parse the slow galvo trace to find frame times relative to analog signals

[frameStarts, incompleteFrames] = parseSlowGalvo(slowGalvo);

%% Load neural data from Suite2P, convert to a more useful format (A and S),
% IGNORE THIS CELL IF USING NEWEST PYTHON VERSION OF SUITE2P!!!!

% if ischar(datFilename)
%   load(datFilename);
% elseif isstruct(datFilename)
%   dat = datFilename;
% else
%   error('Second argument must be data path or dat struct');
% end
% 
% [data.A, S, data.neuropilCoeffs, data.kernels, recon] = suite2PToAS(dat, showPlots);
% data.im = dat.mimg(:, :, dat.maxmap);

% data.neuropilCoeffs, data.kernels will no be used later in the script
% data.A will only be used to count its number of elements
% recon will be used in the next cell, which can be skipped

%% Produce dF/F
% IGNORE THIS CELL IF USING THE LATEST VERSION OF SUITE2P

% dFOF = konnerthDeltaFOverF(recon')';


%% assign spks to S and DS to xoff and yoff

S = spks((iscell(:,1)==1),:);
DS = [ops.xoff' ops.yoff']; % This step is necessary due to suite2p upgrade

%% Check frame count

% MScan usually drops partial frames
nExtrasFrameStarts = length(frameStarts) - size(S, 2);
if length(incompleteFrames) >= nExtrasFrameStarts
  fprintf('Discarding %d incomplete frame(s) from slow galvo trace\n', length(incompleteFrames));
  frameStarts(incompleteFrames(end-nExtrasFrameStarts+1:end)) = [];
end

frameDiff = abs(length(frameStarts) - size(S,2));
if length(frameStarts) > size(S, 2)
  warning(['Parsed galvo has more frames as Suite2P output! Dropped ' num2str(frameDiff) ' to make them equal.']);
  frameStarts = frameStarts(1:size(S,2));
elseif length(frameStarts) < size(S, 2)
  warning(['Parsed galvo has less frames as Suite2P output! Dropped ' num2str(frameDiff) ' to make them equal.']);
  S = S(:, 1:length(frameStarts));
end


%% Find alignment frames

nTrials = length(stimOnSamples);

fr = 1;
stimFrames = NaN(1, nTrials);

for tr = 1:nTrials
    while fr < length(frameStarts) && frameStarts(fr+1) < stimOnSamples(tr)
        fr = fr + 1;
    end
    
    if fr == length(frameStarts)
        error('Ran out of frames');
    end
    
    if abs(stimOnSamples(tr) - frameStarts(fr)) <= abs(stimOnSamples(tr) - frameStarts(fr + 1))
        stimFrames(tr) = fr;
    else
        stimFrames(tr) = fr + 1;
    end
end


%% Pack trial aligned 2P data

preFrames = ceil(1000 * preTime / msPerFrame);
postFrames = ceil(1000 * postTime / msPerFrame);

framesPerTrial = preFrames + postFrames;
data.trialStimFrame = preFrames;
data.stimFramesOrig = stimFrames;

% data.neural = zeros(size(data.A, 2), framesPerTrial, nTrials); %change data.A
data.neural = zeros(sum(iscell(:,1)), framesPerTrial, nTrials); %changed dsize(data.A, 2) to sum(iscell(:,1))
data.neuralTimes = -preFrames + 1 : postFrames;
data.neuralTimesMs = data.neuralTimes * msPerFrame;

% data.dFOF = zeros(size(data.A, 2), framesPerTrial, nTrials); % change data.A
% data.dFOF = zeros(sum(iscell(:,1)), framesPerTrial, nTrials); %changed dsize(data.A, 2) to sum(iscell(:,1))
data.DS = zeros(2, framesPerTrial, nTrials);

for tr = 1:nTrials
    data.neural(:, :, tr) = S(:, stimFrames(tr) + data.neuralTimes);
    data.DS(:, :, tr) = DS(stimFrames(tr) + data.neuralTimes, :)';  % modified line
%   data.dFOF(:, :, tr) = dFOF(:, stimFrames(tr) + data.neuralTimes); 
%   data.DS(:, :, tr) = dat.ops.DS(stimFrames(tr) + data.neuralTimes, :)';
end


%% Pack trial aligned analog data

preSamples = ceil(data.analogFreq * preTime);
postSamples = ceil(data.analogFreq * postTime);

samplesPerTrial = preSamples + postSamples;
data.trialStimSample = preSamples;
data.stimSamplesOrig = NaN(1, length(stimFrames));

data.analog = zeros(length(analogChannels), samplesPerTrial, nTrials);
data.analogTimes = -preSamples + 1 : postSamples;

for tr = 1:nTrials
  stimOnSample = frameStarts(stimFrames(tr));
  data.stimSamplesOrig(tr) = stimOnSample;
  data.analog(:, :, tr) = volt(analogChannels, stimOnSample + data.analogTimes);
end 

%% Save the indices of red cells
data.idx_redcell = find(redcell(find(iscell(:,1)),1));
data.idx_notredcell = find(~redcell(find(iscell(:,1)),1));

%% Save data in data.mat
[baseDir binFileName binExt]= fileparts(binFilename);
save([baseDir filesep 'suite2p' filesep 'plane0' filesep 'data.mat'],'data');
