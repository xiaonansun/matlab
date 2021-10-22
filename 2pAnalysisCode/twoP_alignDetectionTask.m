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


%% Load analog data (MScan .bin file)

tic
if numel(binFilename) > 1
    disp(['There are ' num2str(numel(binFilename)) ' binary MScan files.']);
    for i = 1:numel(binFilename)
        [volt{i}, data.analogFreq{i}] = readMOMAnalog(binFilename{i});
        disp(['MScan analog file (.bin) loaded in ' num2str(toc) ' seconds.']);
        
        trialCodes{i} = volt{i}(trialCodeCh, :);
        slowGalvo{i} = volt{i}(slowGalvoCh, :);
        trialStop{i} = volt{i}(trialOffCh, :);
        stimOn{i} = volt{i}(stimOnCh, :);
    end
else
    [volt, data.analogFreq] = readMOMAnalog(binFilename{1});
    disp(['MScan analog file (.bin) loaded in ' num2str(toc) ' seconds.']);
    
    trialCodes = volt(trialCodeCh, :);
    slowGalvo = volt(slowGalvoCh, :);
    trialStop = volt(trialOffCh, :);
    stimOn = volt(stimOnCh, :);
end

%% Find alignment values, decode trial codes
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        stimOn{i}(stimOn{i} <= voltageThresh) = 0;
        stimOn{i}(stimOn{i} > voltageThresh) = 1;
        stimOnSamples{i} = find(diff(stimOn{i}) == 1) + 1; % stimOnSamples is the stimulus onset index of the analog voltage trace
        
        [trialNumbers{i}, codeStarts{i}] = segmentVoltageAndReadBarcodes(trialCodes{i}, shortInt, longInt);
        
        trialStop{i}(trialStop{i} < voltageThresh) = 0;
        trialStop{i}(trialStop{i} > voltageThresh) = 1;
        trialStopSamples{i} = find(diff(trialStop{i}) == 1) + 1;
        
        if i == 1
            idxLastTrialNumber = max(find(trialNumbers{i} >0));
            stimOnSamples{i} = stimOnSamples{i}(1:end-1);
            trialNumbers{i} = trialNumbers{i}(1:idxLastTrialNumber-1);
            codeStarts{i} = codeStarts{i}(1:idxLastTrialNumber -1);
            trialStopSamples{i} = trialStopSamples{i}(1:end-1);
        end
    end
else
stimOn(stimOn <= voltageThresh) = 0;
stimOn(stimOn > voltageThresh) = 1;
stimOnSamples = find(diff(stimOn) == 1) + 1;

[trialNumbers, codeStarts] = segmentVoltageAndReadBarcodes(trialCodes, shortInt, longInt);

data.trialCodesRaw = trialNumbers;
data.trialCodes = trialNumbers(trialNumbers>0);

trialStop(trialStop < voltageThresh) = 0;
trialStop(trialStop > voltageThresh) = 1;
trialStopSamples = find(diff(trialStop) == 1) + 1;
end


%% check if trialNumbers are monotonic. If not if they could belong to multiple sessions.
% tNCell = trialNumbers; % Assign trialNumbers to another variable, since trialNumbers will be changed
% if iscell(trialNumbers)
%     trialNumbers = horzcat(trialNumbers{:});
% end

if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        trials{i} = trialNumbers{i}(trialNumbers{i} > -1);
        nrChange{i} = find(diff(trials{i}) ~= 1);
        
        if ~isempty(nrChange{i})
            if nrChange{i} > 100 || length(nrChange{i}) > 1
                warning('Trial numbers are not monotonic! Increment subsequent trialnumbers to recover full session!')
                for iChange = 1 : length(nrChange{i})
                    trials{i}(nrChange{i} + 1 : end) = trials{i}(nrChange{i} + 1 : end) + trials{i}(nrChange{i});
                end
                trialNumbers{i}(trialNumbers{i} > -1) = trials{i};
            else
                if length(trials{i}) == length(stimOnSamples{i})
                    warning('Trial numbers are not monotonic after less than 100 trials. Rejecting those trials!');
                    stimOnSamples{i}(1:nrChange{i}) = [];
                else
                    error('Trial numbers are not monotonic and dont match stimulus onset times');
                end
            end
        end
    end
else
    trials = trialNumbers(trialNumbers > -1);
    nrChange = find(diff(trials) ~= 1);
    
    if ~isempty(nrChange)
        %         if nrChange > 100 || length(nrChange) > 1
        if nrChange > 10 || length(nrChange) > 1
            warning('Trial numbers are not monotonic! Increment subsequent trialnumbers to recover full session!')
            for iChange = 1 : length(nrChange)
                trials(nrChange + 1 : end) = trials(nrChange + 1 : end) + trials(nrChange);
            end
            trialNumbers(trialNumbers > -1) = trials;
        else
            if length(trials) == length(stimOnSamples)
                warning('Trial numbers are not monotonic after less than 100 trials. Rejecting those trials!');
                stimOnSamples(1:nrChange) = [];
            else
                error('Trial numbers are not monotonic and dont match stimulus onset times');
            end
        end
    end
end
%% Now, find the starts, stops, and trial codes corresponding to each stimOn
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        [starts{i}, stops{i}, data.trialNumbers{i}] = alignStartsStopsNumbers(stimOnSamples{i}, codeStarts{i}, trialNumbers{i}, trialStopSamples{i});
        if any(isnan(starts{i}))
            stops{i}(isnan(starts{i})) = [];
            data.trialNumbers{i}(isnan(starts{i})) = [];
            starts{i}(isnan(starts{i})) = [];
        end
    end
else
    [starts, stops, data.trialNumbers] = alignStartsStopsNumbers(stimOnSamples, codeStarts, trialNumbers, trialStopSamples);
    if any(isnan(starts))
        stops(isnan(starts)) = [];
        data.trialNumbers(isnan(starts)) = [];
        starts(isnan(starts)) = [];
    end
end
%% Parse the slow galvo trace to find frame times relative to analog signals
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        [frameStarts{i}, incompleteFrames{i}] = parseSlowGalvo(slowGalvo{i});
    end
else
    [frameStarts, incompleteFrames] = parseSlowGalvo(slowGalvo); % frameStarts: length is the total number of frames and values are the indices of the beginning of each frame based on the analog voltage trace
end



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
if numel(binFilename) > 1
    idxLastFiles = find(ops.frames_per_file ~= mode(ops.frames_per_file));
    idxFrame = 1;
    for i = 1:numel(binFilename)
        spksCell{i} = spks(:,idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))));
        S{i} = spksCell{i}((iscell(:,1)==1),:);
        DS{i} = [ops.xoff(:,idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))))' ops.yoff(:,idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))))']; % This step is necessary due to suite2p upgrade
        idxFrame = sum(ops.frames_per_file(1:idxLastFiles(i)))+1;
    end
else
    S = spks((iscell(:,1)==1),:);
    DS = [ops.xoff' ops.yoff']; % This step is necessary due to suite2p upgrade 
end

%% Check frame count

% MScan usually drops partial frames
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        nExtrasFrameStarts{i} = length(frameStarts{i}) - size(S, 2);
        if length(incompleteFrames{i}) >= nExtrasFrameStarts{i}
            fprintf('Discarding %d incomplete frame(s) from slow galvo trace\n', length(incompleteFrames{i}));
            frameStarts{i}(incompleteFrames{i}(end-nExtrasFrameStarts{i}+1:end)) = [];
        end
        
        frameDiff{i} = abs(length(frameStarts{i}) - size(S,2));
        if length(frameStarts{i}) > size(S{i}, 2)
            warning(['Parsed galvo has more frames as Suite2P output! Dropped ' num2str(frameDiff{i}) ' to make them equal.']);
            frameStarts{i} = frameStarts{i}(1:size(S{i},2));
        elseif length(frameStarts{i}) < size(S{i}, 2)
            warning(['Parsed galvo has less frames as Suite2P output! Dropped ' num2str(frameDiff{i}) ' to make them equal.']);
            S{i} = S{i}(:, 1:length(frameStarts{i}));
        end
    end
else
    nExtrasFrameStarts = length(frameStarts) - size(S, 2);
    if length(incompleteFrames) >= nExtrasFrameStarts
        fprintf('Discarding %d incomplete frame(s) from slow galvo trace\n', length(incompleteFrames));
        frameStarts(incompleteFrames(end-nExtrasFrameStarts+1:end)) = [];
    end
    
    frameDiff = abs(length(frameStarts) - size(S,2));
    if length(frameStarts) > size(S, 2)
        warning(['Parsed galvo has more frames than Suite2P output! Dropped ' num2str(frameDiff) ' frames from parsed galvo to make them equal.']);
        frameStarts = frameStarts(1:size(S,2));
    elseif length(frameStarts) < size(S, 2)
        warning(['Parsed galvo has fewer frames than Suite2P output! Dropped ' num2str(frameDiff) ' frames from Suite2P output to make them equal.']);
        S = S(:, 1:length(frameStarts));
    end
end
%% Find alignment frames
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        nTrials{i} = length(stimOnSamples{i});
        fr = 1;
        stimFrames{i} = NaN(1, nTrials{i});
        for tr = 1:nTrials{i}
            while fr < length(frameStarts{i}) && frameStarts{i}(fr+1) < stimOnSamples{i}(tr)
                fr = fr + 1;
            end
            
            if fr == length(frameStarts{i})
                error('Ran out of frames');
            end
            
            if abs(stimOnSamples{i}(tr) - frameStarts{i}(fr)) <= abs(stimOnSamples{i}(tr) - frameStarts{i}(fr + 1))
                stimFrames{i}(tr) = fr;
            else
                stimFrames{i}(tr) = fr + 1;
            end
        end
    end
else
    nTrials = length(stimOnSamples);
    
    fr = 1;
    stimFrames = NaN(1, nTrials);
    
    for tr = 1:nTrials
        while fr < length(frameStarts) && frameStarts(fr+1) < stimOnSamples(tr)
            fr = fr + 1;
        end
        
        if fr == length(frameStarts)
            error('Ran out of frames: check the number of frames in 2P tiff files as well as in spks');
        end
        
        if abs(stimOnSamples(tr) - frameStarts(fr)) <= abs(stimOnSamples(tr) - frameStarts(fr + 1))
            stimFrames(tr) = fr;
        else
            stimFrames(tr) = fr + 1;
        end
    end
end
%% Pack trial aligned 2P data
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        preFrames = ceil(1000 * preTime / msPerFrame);
        postFrames = ceil(1000 * postTime / msPerFrame);
        
        framesPerTrial = preFrames + postFrames;
        data.trialStimFrame{i} = preFrames;
        data.stimFramesOrig{i} = stimFrames{i};
        
        data.neural{i} = zeros(sum(iscell(:,1)), framesPerTrial, nTrials{i}); %changed dsize(data.A, 2) to sum(iscell(:,1))
        data.neuralTimes{i} = -preFrames + 1 : postFrames;
        data.neuralTimesMs{i} = data.neuralTimes{i} * msPerFrame;
        
        data.DS{i} = zeros(2, framesPerTrial, nTrials{i});
        
        for tr = 1:nTrials{i}
            data.neural{i}(:, :, tr) = S{i}(:, stimFrames{i}(tr) + data.neuralTimes{i});  %% ISSUE HERE: the last trial before the first recording was terminated is incomplete, leading to a "Index in position 2 exceeds array bounds" error. Will need to remove this trial early in the code.
            data.DS{i}(:, :, tr) = DS{i}(stimFrames{i}(tr) + data.neuralTimes{i}, :)';  % modified line
        end
    end
else
    preFrames = ceil(1000 * preTime / msPerFrame);
    postFrames = ceil(1000 * postTime / msPerFrame);
    
    framesPerTrial = preFrames + postFrames;
    data.trialStimFrame = preFrames;
    data.stimFramesOrig = stimFrames;
    
    data.neural = zeros(sum(iscell(:,1)), framesPerTrial, nTrials); %changed dsize(data.A, 2) to sum(iscell(:,1))
    data.neuralTimes = -preFrames + 1 : postFrames;
    data.neuralTimesMs = data.neuralTimes * msPerFrame;
    
    data.DS = zeros(2, framesPerTrial, nTrials);
    
    for tr = 1:nTrials
        try
        data.neural(:, :, tr) = S(:, stimFrames(tr) + data.neuralTimes);
        data.DS(:, :, tr) = DS(stimFrames(tr) + data.neuralTimes, :)';  % modified line
        catch ME
            disp(ME.message);
        end
    end
end
%% Pack trial aligned analog data
if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        preSamples{i} = ceil(data.analogFreq{i} * preTime);
        postSamples{i} = ceil(data.analogFreq{i} * postTime);
        
        samplesPerTrial{i} = preSamples{i} + postSamples{i};
        data.trialStimSample{i} = preSamples{i};
        data.stimSamplesOrig{i} = NaN(1, length(stimFrames{i}));
        
        data.analog{i} = zeros(length(analogChannels), samplesPerTrial{i}, nTrials{i});
        data.analogTimes{i} = -preSamples{i} + 1 : postSamples{i};
        
        for tr = 1:nTrials{i}
            stimOnSample{i} = frameStarts{i}(stimFrames{i}(tr));
            data.stimSamplesOrig{i}(tr) = stimOnSample{i};
            data.analog{i}(:, :, tr) = volt{i}(analogChannels, stimOnSample{i} + data.analogTimes{i});
        end
    end
else
    preSamples = ceil(data.analogFreq * preTime);
    postSamples = ceil(data.analogFreq * postTime);
    
    samplesPerTrial = preSamples + postSamples;
    data.trialStimSample = preSamples;
    data.stimSamplesOrig = NaN(1, length(stimFrames));
    
    data.analog = zeros(length(analogChannels), samplesPerTrial, nTrials);
    data.analogTimes = -preSamples + 1 : postSamples;
    
    for tr = 1:nTrials
        try
        stimOnSample = frameStarts(stimFrames(tr));
        data.stimSamplesOrig(tr) = stimOnSample;
        data.analog(:, :, tr) = volt(analogChannels, stimOnSample + data.analogTimes);
        catch ME
            disp(ME.message);
        end
    end
end
%% Save the indices of red cells
data.idx_redcell = find(redcell(find(iscell(:,1)),1));
data.idx_notredcell = find(~redcell(find(iscell(:,1)),1));

%% Save data in data.mat
[baseDir,~,~]= fileparts(binFilename{1});
save([baseDir filesep 'suite2p' filesep 'plane0' filesep 'data.mat'],'data');
