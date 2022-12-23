function data = twoP_alignDetectionTask(ops, spks, isneuron, redcell, binFilename,showPlots)
% data = twoP_alignDetectionTask(datFilename, binFilename [, showPlots])
% data = twoP_alignDetectionTask(dat, binFilename [, showPlots])
% data = twoP_alignDetectionTask(ops, spks, iscell, redcell, binFilename, showPlots)

% ----- 2021 Revision -----=
% This code now accomodates for sessions with multiple MScan binary files,
% which indicates that 2p acquisition was stopped and then re-started,
% usually due to an MScan computer crash. 

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

%% commented variables

% animal = 'CSP27'; session = '20200319a';
% npy = twoP_importSuite2pNPY(animal,session); % use this when testing script
% ops = npy.ops; spks = npy.spks; iscell=npy.iscell; redcell=npy.redcell;
% binFilename= npy.bin_MScan_filepath;

%% Parameters

% Sampling

% This value is shown (as frame duration) in the File Properties button of
% File Manager window of Mview, presumably more accurate than 1000/30.9 Hz
msPerFrame =  32.363833;

data.msPerFrame = msPerFrame;

% Time window (in seconds)
preTime = 3; % time window prior to the stimulus included in PSTH 
postTime = 4; % time window after the stimulus included in PSTH 

% TTL interpretation
voltageThresh = 1.5;  % V
% Simon used longer intervals for barcodes than Matt Kaufman's defaults
shortInt = 2.5;
longInt = 5.5;

% assign analog channels. These are the channels recorded in the binary
% files from MScan
trialCodeCh = 1; % trial code (aka trial ID numer) in analog TTL binary form
slowGalvoCh = 2; % analog signal synchronized with the 2p galvo for tracking frames
trialOffCh = 3; % Trial end TTL sent to MSscan data acquisition from Bpod
stimOnCh = 4; % Stimulus onset TTL sent to MSscan data acquisition from Bpod

analogChannels = 1:5;  % channels to save
data.analogIdentities = {'trialCode', 'slowGalvo', 'trialOff', 'stimOn', 'piezo'};

data.animal = ops.animal;
data.session = ops.session;

%% Computes SDU traces
overwriteSDU = false;
sdu_trace = twoP_computeZScoreFromSpks(ops.animal,ops.session,overwriteSDU);

%% Optional arguments

if ~exist('showPlots', 'var')
  showPlots = 0;
end

%% Load analog data (MScan .bin file)

tic
if numel(binFilename) > 1 %if there are more than one MSCan binary files, then the imaging session was interrupted by an MSCan crash
    disp(['There are ' num2str(numel(binFilename)) ' binary MScan files.']);
    for i = 1:numel(binFilename)
        [volt{i}, data.analogFreq{i}] = readMOMAnalog(binFilename{i}); % IMPORTANT! This function reads MScan binary files
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
    trialStop = volt(trialOffCh, :); % Analog voltage of the stop TTL signal. Maximum value ~ 4.8, minimum value ~ -0.0037.  sum(trialStop < 3 & trialStop > 0.1) = 0
    stimOn = volt(stimOnCh, :); % analog voltage of the stimulus-onset TTL. sum(stimOn < 2 & stimOn > 0.1) = 0
end

%% Find alignment values, decode trial codes

% This function (segmentVoltageAndReadBarcodes.m) reads trial codes from
% the MScan binary

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
stimOnSamples = find(diff(stimOn) == 1) + 1; %Analog voltage trace index of stimulus onset

[trialNumbers, codeStarts] = segmentVoltageAndReadBarcodes(trialCodes, shortInt, longInt); 
% trialNumber is the trial IDs decoded from the trialCodes analog voltage trace
% codeStarts is the voltage trace index of the arrival of the trial code

data.trialCodesRaw = trialNumbers;
data.trialCodes = trialNumbers(trialNumbers>0);

trialStop(trialStop < voltageThresh) = 0;
trialStop(trialStop > voltageThresh) = 1;
trialStopSamples = find(diff(trialStop) == 1) + 1; % Analog index of trial stop signal, numel(trialStopSamples) yields number of stop signals
end


%% check if trialNumbers are monotonic. If not if they could belong to multiple sessions.
% For example: if the trialNumbers drops down to 1, then it is possible
% Bpod has crashed/reset, hence starting a new behavior session.
% If trialNumbers jumps, especially by more than 10, then it is possible
% the MScan computer crashed and was reset while the bpod computer kept on
% running

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
    trials = trialNumbers(trialNumbers > -1); % This is the vector of trial IDs without bad TTL trial codes
    nrChange = find(diff(trials) ~= 1); % Find trials when the trialCodes sent by Bpod jumps by more than 1.
    
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
% uses this script: alignStartsStopsNumbers.m

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
% uses this script: parseSlowGalvo.m
% Outputs:
% frameStarts: (1) length is the total number of frames counted from the galvo
% voltage trace, and (2) values are the indices of the beginning of each 
%frame based on the analog voltage trace

if numel(binFilename) > 1
    for i = 1:numel(binFilename)
        [frameStarts{i}, incompleteFrames{i}] = parseSlowGalvo(slowGalvo{i});
    end
else
    [frameStarts, incompleteFrames] = parseSlowGalvo(slowGalvo); 
end

% Can compare length(framStarts) with length(spks). If the two are equal,
% then the suite2p-processed imaging data contains the same number of
% frames as those counted from the galvo analog signal. This suggests
% consistency in data processing. If the two are different by 1 frame, then
% the last frame can be ignored. The extra "frame" is usually on the end of
% the slow galvo voltage trace.

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
% S is the deconvolved spkes of ROIs identified as real cells

if numel(binFilename) > 1
    idxLastFiles = find(ops.frames_per_file ~= mode(ops.frames_per_file));
    idxFrame = 1;
    for i = 1:numel(binFilename)
%         spksCell{i} = spks((isneuron(:,1)==1),idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))));
        S{i} = spks((isneuron(:,1)==1),idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))));
%         S{i} = spksCell{i}((isneuron(:,1)==1),:);
        SDU{i} = sdu_trace((isneuron(:,1)==1),idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))));
        DS{i} = [ops.xoff(:,idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))))' ops.yoff(:,idxFrame:sum(ops.frames_per_file(1:idxLastFiles(i))))']; % This step is necessary due to suite2p upgrade
        idxFrame = sum(ops.frames_per_file(1:idxLastFiles(i)))+1;
    end
else
    S = spks((isneuron(:,1)==1),:);
    SDU = sdu_trace((isneuron(:,1)==1),:);
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
    nExtrasFrameStarts = length(frameStarts) - size(S, 2); % The number of extra frames: difference between the number of frames counted from the galvo voltage trace and the number of frames counted from deconvolved the spike trace
    if length(incompleteFrames) >= nExtrasFrameStarts
        fprintf('Discarding %d incomplete frame(s) from slow galvo trace\n', length(incompleteFrames));
        frameStarts(incompleteFrames(end-nExtrasFrameStarts+1:end)) = []; % Removes the extra frames from the slow galvo voltage trace
    end
    
    % The next few lines assumes that if there are extra frames recorded on
    % either the slow galvo end or from the imaging end, those frames recorded 
    % at the END of the session and can be safely removed to match the
    % frame count of the other, whichever is less.
    
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
    
    for tr = 1:nTrials % loop through all trials with a TTL stimulus-on signal
        while fr < length(frameStarts) && frameStarts(fr+1) < stimOnSamples(tr)
            fr = fr + 1;
        end
        % In this while loop, the first condition tests for the current
        % frame index to stay within the total frame count of the session.
        % the second condition tests for the MScan voltage trace index and
        % stops when fr arrives at a stimulus on TTL signal.
        
        if fr == length(frameStarts)
            error('Ran out of frames: check the number of frames in 2P tiff files as well as in spks');
        end
        % This condition will yield an error if the current frame count
        % exceeds the total number of images collected.
        
        if abs(stimOnSamples(tr) - frameStarts(fr)) <= abs(stimOnSamples(tr) - frameStarts(fr + 1))
            stimFrames(tr) = fr;
        else
            stimFrames(tr) = fr + 1;
        end
        % This condition compares the difference between the stim-on
        % voltage index of the Nth trial and the current frame's voltage
        % index, versus the difference between the stim-on
        % voltage index of the Nth trial and the next frame's voltage
        % index. If the condition is equal or less, then stim-on occured
        % during the current frame. If the condition is greater, then
        % stim-on occured during the subsequent frame.
        
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
        
        data.neural{i} = zeros(sum(isneuron(:,1)), framesPerTrial, nTrials{i}); %changed dsize(data.A, 2) to sum(iscell(:,1))
        data.sdu{i} = zeros(sum(isneuron(:,1)), framesPerTrial, nTrials{i}); %changed dsize(data.A, 2) to sum(iscell(:,1))
        data.neuralTimes{i} = -preFrames + 1 : postFrames;
        data.neuralTimesMs{i} = data.neuralTimes{i} * msPerFrame;
        
        data.DS{i} = zeros(2, framesPerTrial, nTrials{i});
        
        for tr = 1:nTrials{i}
            data.neural{i}(:, :, tr) = S{i}(:, stimFrames{i}(tr) + data.neuralTimes{i});  %% ISSUE HERE: the last trial before the first recording was terminated is incomplete, leading to a "Index in position 2 exceeds array bounds" error. Will need to remove this trial early in the code.
            data.sdu{i}(:, :, tr) = SDU{i}(:, stimFrames{i}(tr) + data.neuralTimes{i});  %% ISSUE HERE: the last trial before the first recording was terminated is incomplete, leading to a "Index in position 2 exceeds array bounds" error. Will need to remove this trial early in the code.
            data.DS{i}(:, :, tr) = DS{i}(stimFrames{i}(tr) + data.neuralTimes{i}, :)';  % modified line
        end
    end
else
    preFrames = ceil(1000 * preTime / msPerFrame);
    postFrames = ceil(1000 * postTime / msPerFrame);
    
    framesPerTrial = preFrames + postFrames;
    data.trialStimFrame = preFrames;
    data.stimFramesOrig = stimFrames;
    data.neuralTimes = -preFrames + 1 : postFrames;
    data.neuralTimesMs = data.neuralTimes * msPerFrame;
    data.DS = zeros(2, framesPerTrial, nTrials);
    data.neural = zeros(sum(isneuron(:,1)), framesPerTrial, nTrials); %changed dsize(data.A, 2) to sum(iscell(:,1))
    data.sdu = zeros(sum(isneuron(:,1)), framesPerTrial, nTrials); %changed dsize(data.A, 2) to sum(iscell(:,1))
    
    for tr = 1:nTrials
        try
        data.neural(:, :, tr) = S(:, stimFrames(tr) + data.neuralTimes);
        data.sdu(:, :, tr) = SDU(:, stimFrames(tr) + data.neuralTimes);
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
data.idx_redcell = find(redcell(find(isneuron(:,1)),1));
data.idx_notredcell = find(~redcell(find(isneuron(:,1)),1));

%% Save data in data.mat
[baseDir,~,~]= fileparts(binFilename{1});
save([baseDir filesep 'suite2p' filesep 'plane0' filesep 'data.mat'],'data');
