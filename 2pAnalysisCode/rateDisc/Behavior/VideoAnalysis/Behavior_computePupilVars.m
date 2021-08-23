function Behavior_computePupilVars(fPath,reload,invertPupil)
% Needs eyeTrace.mj2 videos for pupil. This code is meant to recompute
% pupil diameter and will also produce a filtered version of pupil diamter
% to aproximate modulation based on ACh and NE influence.

% fPath = 'U:\space_managed_data\BpodImager\Animals\mSM43\SpatialDisc\21-Nov-2017\BehaviorVideo\';  %example server data path

if ~exist('reload', 'var')
    reload = true; %recompute pupil diamter
end

if ~exist('invertPupil', 'var')
    invertPupil = false; %assume dark pupil by default
end

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end

if exist([fPath 'bhvOpts.mat'],'file') %check bhvOpts and load if present
    load([fPath 'bhvOpts']);
else
    bhvOpts.maxFrameCnt = 500; %max number of frames per trial. Earlier frames will be removed from analysis.
end

bhvOpts.lowCut = 0.03; %lower cutoff for low-pass filter (for ACh)
bhvOpts.highCut = 10; %higher cutoff for band-pass filter (for NE)
bhvOpts.filterOrder = 50; %filter order for FIR filter
bhvOpts.preStimDur = 3; %baseline duration for trial-aligned face data
highTresh = 0.25;
lowThresh = 0;

%% check files and load SVD data
eyeFiles = dir([fPath '*eyeTrace*.mj2']);

%% loop through trials and perform analysis
tic;
for iTrials = 1:length(eyeFiles)
    
    load([fPath 'faceVars_' int2str(iTrials) '.mat']); %load last dataset to get cTimes
    frameFile = dir([fPath '*frameTimes_'  num2str(iTrials,'%04i') '*_1.mat']); %load file with absolute timestamps
    load([fPath frameFile.name]);
    
    if length(frameTimes) > bhvOpts.maxFrameCnt
        frameTimes = frameTimes(end-bhvOpts.maxFrameCnt+1:end);
    end
    
    if length(cTimes) ~= length(frameTimes) %check if relative and absolute timestamps have the same length. If not there might be a mixup with the files that are used.
        error(['Different framecount for relative and absolute timestamps. Current file: ' frameFile.name])
    end
    
    if reload || isempty(eyeVars)
        eyeTrace = squeeze(importdata([fPath 'eyeTrace_' int2str(iTrials) '.mj2'])); %load eye data for current trial
        eyeTrace = mat2gray(eyeTrace);
        if invertPupil
            %pupil is bright instead of dark. this happens during 2p imaging.
            temp = mat2gray(abs(fft(eyeTrace,[],3))); %use fft to get a better idea what is pupil and what is IR light reflection
            mask = mat2gray(max(temp(:,:,round(size(eyeTrace,3)/2 - 15):round(size(eyeTrace,3)/2 + 15)),[],3)); %find oscillation around half the collected framerate. should be due to 2p laser.
            mask = mask < 0.1; %reject areas that barely oscillate
            mask = mask | temp(:,:,1) > 0.75; %reject area with high 0-response
            eyeTrace = 1 - eyeTrace;
            eyeTrace = reshape(eyeTrace,[],size(eyeTrace,3));
            eyeTrace(mask(:),:) = 1;
            eyeTrace = reshape(eyeTrace,size(mask,1),size(mask,2),[]);
            highTresh = 0.6;
            lowThresh = 0.4;
        end
        
        eyeVars.orient = zeros(size(eyeTrace,3),1);
        eyeVars.center = zeros(size(eyeTrace,3),2);
        eyeVars.axes = zeros(size(eyeTrace,3),2);
        eyeVars.solidity = zeros(size(eyeTrace,3),1);
        
        for iFrames = 1:size(eyeTrace,3)
            temp = Behavior_EyeCheck(eyeTrace(:,:,iFrames),highTresh,lowThresh,false); %get pupil data
            eyeVars.center(iFrames,:) = temp.center; %circle center
            eyeVars.axes(iFrames,:) = temp.axes; %circle axis
            eyeVars.orient(iFrames) = temp.orient; %circle orientation
            eyeVars.solidity(iFrames) = temp.solidity; %circle solidity
        end
        
        % find bad pupil estimates and try to repair.
        trace = mean(eyeVars.axes,2);
        if sum(trace == 0) ~= length(trace)
            trace(trace == 0) = NaN;
            trace = fillgaps(trace, 100);
            smoothTrace = smooth(trace,10,'rlowess'); %use smooth pupil trace as reference for outliers
        else
            smoothTrace = trace;
        end
        smoothDiff = abs(mean(eyeVars.axes,2)-smoothTrace); %difference from smoothed pupil trace
        smoothStd = std(smoothDiff);
        smoothDiff = smoothDiff ./ smoothStd; %put to STDs
        idx = find([0;diff(zscore(mean(eyeVars.axes,2)))]<-.5 | [0;diff(zscore(mean(eyeVars.axes,2)))]>.5 | smoothDiff > 0.5); %potentially bad results. Try to do better by finding result that is closest to average.

        for iFrames = 1:length(idx)
            clear temp tDiff
            Cnt = 0;
            for x = [-6 -2 2 6] %try eroding or fusing patches to get better pupil estimate
                Cnt = Cnt + 1;
                temp{Cnt} = Behavior_EyeCheck(eyeTrace(:,:,idx(iFrames)),highTresh,lowThresh,false,true,x); %get pupil data. Force individual areas to fuse together.
                tDiff(Cnt) = abs(smoothTrace(idx(iFrames)) - mean(temp{Cnt}.axes,2)) ./ smoothStd; %find modification that gets estimate closer to smoothed trace
            end
            [~,minDiff] = min(tDiff);
            if tDiff(minDiff) < abs(trace(idx(iFrames))-smoothTrace(idx(iFrames)))
                eyeVars.center(idx(iFrames),:) = temp{minDiff}.center; %circle center
                eyeVars.axes(idx(iFrames),:) = temp{minDiff}.axes; %circle axis
                eyeVars.orient(idx(iFrames)) = temp{minDiff}.orient; %circle orientation
                eyeVars.solidity(idx(iFrames)) = temp{minDiff}.solidity; %circle solidity
            end
        end
        save([fPath 'faceVars_' int2str(iTrials) '.mat'],'eyeVars','-append');
    end
    
    pupil{iTrials} = mean(eyeVars.axes,2); %keep pupil trace
    snout{iTrials} = snoutMotion'; %keep snout trace
    sniff{iTrials} = noseMotion'; %keep sniff trace
    face{iTrials} = faceMotion'; %keep face trace
    
    pupilTime{iTrials} = frameTimes * 86400; %keep absolute time in seconds
    trialStart(iTrials) = pupilTime{iTrials}(1);
    clear snoutMotion noseMotion eyeVars frameTimes
    
    % load data from body cam
    bodyM{iTrials} = []; bTime{iTrials} = [];
    try
        load([fPath 'bodyVars_' int2str(iTrials) '.mat']); %load last dataset to get cTimes
        frameFile = dir([fPath '*frameTimes_'  num2str(iTrials,'%04i') '*_2.mat']); %load file with absolute timestamps
        load([fPath frameFile.name]);
        
        if length(frameTimes) > bhvOpts.maxFrameCnt
            frameTimes = frameTimes(end-bhvOpts.maxFrameCnt+1:end);
        end
        if length(cTimes) ~= length(frameTimes) %check if relative and absolute timestamps have the same length. If not there might be a mixup with the files that are used.
            error(['Different framecount for relative and absolute timestamps. Current file: ' frameFile.name])
        end
        bodyM{iTrials} = bodyMotion'; %keep body trace
        bTime{iTrials} = frameTimes * 86400; %keep absolute time in seconds
    end
    
    if iTrials == 1
        sRate = mean(diff(pupilTime{iTrials})); %time between video frames
    end
    
    if iTrials > 1
        timeDiff(iTrials) = pupilTime{iTrials}(1) - trialEnd(iTrials - 1); %time difference between end of last and start of current trial
        extraFrames = floor(timeDiff(iTrials) / sRate); %number of additional frames that are required
        
        pupilTime{iTrials} = [pupilTime{iTrials - 1}(end) + (sRate:sRate:sRate*extraFrames)'; pupilTime{iTrials}]; %fill up time vector with extra frames
        pupil{iTrials} = [NaN(extraFrames,1); pupil{iTrials}]; %fill up pupil vector with extra frames
        snout{iTrials} = [NaN(extraFrames,1); snout{iTrials}]; %fill up snout vector with extra frames
        sniff{iTrials} = [NaN(extraFrames,1); sniff{iTrials}]; %fill up sniff vector with extra frames
        face{iTrials} = [NaN(extraFrames,1); face{iTrials}]; %fill up face vector with extra frames
    end
    trialEnd(iTrials) = pupilTime{iTrials}(end);
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,length(eyeFiles));
        toc
    end
end

%% combine pupil diameter into larger vector
pupil = cat(1,pupil{:});
snout = cat(1,snout{:});
sniff = cat(1,sniff{:});
face = cat(1,face{:});
pupilTime = cat(1,pupilTime{:});

% remove outliers
temp = (pupil - nanmean(pupil)) ./ nanstd(pupil);
idx = [0;diff(temp) > 0.25] | [0;diff(temp) < -0.25] | temp < -4 | temp > 4; %index for blinks / false reads
pupil(idx) = NaN;
pupilFill = fillgaps(pupil, 300); %fill NaN gaps between trials
pupilFill = [pupilFill;pupilFill(end-round(length(pupil) * 0.1):end)]; %add some padding to the end to avoid filter artefacts
pupilFill = zscore(pupilFill);

% high-pass butter filter for pupil data. this is to isolate adrenergic component
fastFilter  = fdesign.bandpass('N,F3dB1,F3dB2', bhvOpts.filterOrder, bhvOpts.lowCut, bhvOpts.highCut, 1/sRate);
fastFilter = design(fastFilter, 'butter');
fastPupil = fastFilter.filter(pupilFill); fastPupil = flipud(fastFilter.filter(flipud(fastPupil)));

% low-pass butter filter for pupil data. this is to isolate cholinergic component
slowFilter  = fdesign.lowpass('N,F3dB', bhvOpts.filterOrder, bhvOpts.lowCut, 1/sRate);
slowFilter = design(slowFilter, 'butter');
slowPupil = slowFilter.filter(pupilFill); slowPupil = flipud(slowFilter.filter(flipud(slowPupil)));

%remove padding
pupilFill(end-round(length(pupil) * 0.1):end) = [];
fastPupil(end-round(length(pupil) * 0.1):end) = [];
slowPupil(end-round(length(pupil) * 0.1):end) = [];

%% reconstruct individual trials and save
orgPupil = cell(1,length(trialStart));
fillPupil = cell(1,length(trialStart));
sPupil = cell(1,length(trialStart));
fPupil = cell(1,length(trialStart));
pTime = cell(1,length(trialStart));
whisker = cell(1,length(trialStart));
nose = cell(1,length(trialStart));
faceM = cell(1,length(trialStart));

for iTrials = 1:length(trialStart)
    idx = pupilTime >= trialStart(iTrials) & pupilTime <= trialEnd(iTrials);
    
    if sum(idx) == 0
        error('Something went wrong with trial indexing.')
    end
    orgPupil{iTrials} = pupil(idx); %original pupil
    fillPupil{iTrials} = pupilFill(idx); %filled pupil
    sPupil{iTrials} = slowPupil(idx); %slow pupil
    fPupil{iTrials} = fastPupil(idx); %fast pupil
    pTime{iTrials} = pupilTime(idx); %timestamps
    whisker{iTrials} = snout(idx); %whisker (snout) motion - this is to ensure that pupil and whisker data have the same timing
    nose{iTrials} = sniff(idx); %nose motion
    faceM{iTrials} = face(idx); %face motion
        
end

save([fPath 'FilteredPupil.mat'],'orgPupil','fillPupil','sPupil','fPupil','pTime','whisker','nose','faceM','bodyM','bTime','pupil','pupilFill','fastPupil','slowPupil','pupilTime');