function Behavior_recomputeEyeTrace(cPath, Animal, Rec, eyeCam)
% Code to convert face variables from compressed movies. Will compute face
% motion, body motion and pupil diameter and save in faceVars.mat.
% This code is meant to be used after converting videos to SVD. 
% Needs eyeTrace.mj2 videos for pupil and V and U components for face variables.

if ~exist('eyeCam', 'var')
    eyeCam = 1; %camera that is viewing the face. Default is cam1.
end

Paradigm = 'SpatialDisc';
fPath = [cPath Animal filesep Paradigm filesep Rec '\BehaviorVideo' filesep]; %server data path

bhvOpts.varCnt = 50; %number of trials used to compute variance map
bhvOpts.eyeFrame = 10; %add some pixels to eyeFrame to ensure pupil is properly captured
bhvOpts.eyeThresh = .25;
bhvOpts.eyeErrode = 2;
bhvOpts.maxFrameCnt = 500; %max number of frames per trial. Earlier frames will be removed from analysis.

%% check files and load SVD data
timeFiles = dir([fPath '*frameTimes*_' int2str(eyeCam) '.mat']); %all files for current cam based on integer before .mat
movieFiles = dir([fPath '*Video*_' int2str(eyeCam) '.mp4']); %all files for current cam based on integer before .mat

bhvFile = dir([fileparts(fPath(1:end-1)) '\' Animal '*.mat']);
load([fileparts(fPath(1:end-1)) '\' bhvFile.name]); %load bpod data

SessionData = Behavior_addNose(fPath, SessionData, eyeCam); %add new indicators
save([fileparts(fPath(1:end-1)) '\' bhvFile.name], 'SessionData'); %save bpod data

eyePos = SessionData.eyePos;
snoutPos = SessionData.snoutPos;
nosePos = SessionData.nosePos;

%% loop through trials and perform analysis
tic;
for iFiles = 1:length(movieFiles)
    %% get raw video and timestamps
    cFile = [fPath movieFiles(iFiles).name];
    rawData = squeeze(importdata(cFile));
    if length(size(rawData)) == 4
        rawData = squeeze(rawData(:, :, 1, :));
    end
    if size(rawData,3) > bhvOpts.maxFrameCnt
        rawData = rawData(:,:,end-bhvOpts.maxFrameCnt+1:end);
    end
    
    load([fPath timeFiles(iFiles).name]) % load frame times
    % check if timestamps are shifted by an hour. Apparently that can happen.
    timeCheck = (SessionData.TrialStartTime(iFiles)*86400) - (frameTimes(1:10) * 86400); %time difference between first acquired frame and onset of current trial
    if any(timeCheck > 3540 & timeCheck < 3660) %timeshift by one hour (+- 10seconds)
        if iFiles == 1
            warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
        end
        frameTimes = frameTimes + (1/24); %add one hour
    elseif any(timeCheck > 30)
        error(['Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds in file: ' timeFiles(iFiles).name])
    end
    
    if size(frameTimes,1) > bhvOpts.maxFrameCnt
        frameTimes = frameTimes(end-bhvOpts.maxFrameCnt+1:end);
    end
    cTimes = ((frameTimes - SessionData.TrialStartTime(iFiles)) * 86400); %camera time stamps, relative to trial onset

    %% isolate eyes, snout and nose
    eyeTrace = rawData((eyePos(2)-eyePos(3)-(bhvOpts.eyeFrame/2)-1) + (1:eyePos(3)*2+bhvOpts.eyeFrame+1), (eyePos(1)-eyePos(3)-(bhvOpts.eyeFrame/2)-1) + (1:eyePos(3)*2+bhvOpts.eyeFrame+1), :);
    eyeTrace = mat2gray(eyeTrace);
    eyeTrace = im2uint8(reshape(eyeTrace,size(eyeTrace,1),size(eyeTrace,2),1,[]));
    v = VideoWriter([fPath 'eyeTrace_' int2str(iFiles) '.mj2'],'Archival'); %write small eye video
    open(v); writeVideo(v,eyeTrace); close(v);
    
    % compute snout motion
    idx1 = (snoutPos(2)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
    idx2 = (snoutPos(1)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
    if any(idx1 > size(rawData,1)) || any(idx2 > size(rawData,2))
        if iFiles == 1
            warning('Snout index is too large. Clipped to fit into image.');
        end
        idx1(idx1 > size(rawData,1)) = [];
        idx2(idx2 > size(rawData,2)) = [];
    end
    snoutMotion = rawData(idx1, idx2, :);
    snoutMotion = diff(snoutMotion,1,3);
    snoutMotion = reshape(snoutMotion,[],size(rawData,3)-1);
    snoutMotion = mean(snoutMotion);
    snoutMotion = [snoutMotion(1) snoutMotion];
    
    %coompute face motion
    faceMotion = diff(rawData,1,3);
    faceMotion = reshape(faceMotion,[],size(rawData,3)-1);
    faceMotion = mean(faceMotion);
    faceMotion = [faceMotion(1) faceMotion];
    
    %compute nose motion
    idx1 = (nosePos(2)-nosePos(3)-1) + (1:nosePos(3)*2+1);
    idx2 = (nosePos(1)-nosePos(3)-1) + (1:nosePos(3)*2+1);
    if any(idx1 > size(rawData,1)) || any(idx2 > size(rawData,2))
        if iFiles == 1
            warning('Nose index is too large. Clipped to fit into image.');
        end
        idx1(idx1 > size(rawData,1)) = [];
        idx2(idx2 > size(rawData,2)) = [];
    end
    noseMotion = rawData(idx1, idx2, :);
    noseMotion = diff(noseMotion,1,3);
    noseMotion = reshape(noseMotion,[],size(rawData,3)-1);
    noseMotion = mean(noseMotion);
    noseMotion = [noseMotion(1) noseMotion];
                                        
    eyeVars = []; %pupil is now computed later
    save([fPath 'faceVars_' int2str(iFiles) '.mat'],'eyeVars','snoutMotion','faceMotion','noseMotion','cTimes');
        
    % give some feedback over progress
    if rem(iFiles,50) == 0
        fprintf(1, 'Current file is %d out of %d\n', iFiles,length(movieFiles));
        toc
    end
end