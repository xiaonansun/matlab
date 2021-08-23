function Behavior_computeBodyVars(cPath, Animal, Rec, bodyCam)
% Code to convert face variables from compressed movies. Will compute face
% motion, body motion and pupil diameter and save in faceVars.mat.
% This code is meant to be used after converting videos to SVD. 
% Needs eyeTrace.mj2 videos for pupil and V and U components for face variables.

if ~exist('bodyCam', 'var')
    bodyCam = 2; %camera that is viewing the body. Default is cam2.
end

Paradigm = 'SpatialDisc';
fPath = [cPath Animal filesep Paradigm filesep Rec '\BehaviorVideo']; %server data path

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end

if exist([fPath 'bhvOpts.mat'],'file') %check bhvOpts and load if present
    load([fPath 'bhvOpts']);
else
    bhvOpts.maxFrameCnt = 500; %max number of frames per trial. Earlier frames will be removed from analysis.
end

%% check files and load raw data
bodyFiles = dir([fPath Animal '*_' num2str(bodyCam) '.mj2']);
timeFiles = dir([fPath '*frameTimes*_' num2str(bodyCam) '.mat']);
bhvFile = dir([fileparts(fPath(1:end-1)) '\' Animal '_' Paradigm '*.mat']);
load([fileparts(fPath(1:end-1)) '\' bhvFile.name]); %load bpod data

%% loop through trials and get body movements
tic;

for iTrials = 1 : min([length(bodyFiles) length(SessionData.TrialStartTime)])
    
    cFile = [fPath bodyFiles(iTrials).name]; %current video file
    rawData = squeeze(importdata(cFile)); %load video data
    load([fPath timeFiles(iTrials).name]); %load time stamps
    
    %check timestamps and number of frames matches up
    if size(rawData,3) ~= length(frameTimes)
        error('Number of frames and frameTimes is unequal. Thats not good.')
    end

    if size(rawData,3) > bhvOpts.maxFrameCnt
        rawData = rawData(:,:,end-bhvOpts.maxFrameCnt+1:end);
        frameTimes = frameTimes(end-bhvOpts.maxFrameCnt+1:end);
    end

    if iTrials == 1
        % check if timestamps are shifted by an hour. Apparently that can happen.
        timeCheck = (SessionData.TrialStartTime(iTrials)*86400) - (frameTimes(1) * 86400); %time difference between first acquired frame and onset of first trial
        if timeCheck > 3590 && timeCheck < 3610 %timeshift by one hour (+- 10seconds)
            if iTrials == 1
                warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
            end
            timeCheck = true;
        elseif timeCheck > 30
            error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
        else
            timeCheck = false;
        end
    end
    if timeCheck
        frameTimes = frameTimes + (1/24); %add one hour
    end
    cTimes = ((frameTimes - SessionData.TrialStartTime(iTrials)) * 86400); %compute timestamps relative to trial onset time
    
    % compute body motion
    bodyMotion = diff(rawData,1,3);
    bodyMotion = reshape(bodyMotion,[],size(rawData,3)-1);
    bodyMotion = mean(bodyMotion);
    bodyMotion = [bodyMotion(1) bodyMotion];
    
    % save data and increase counter
    save([fPath 'bodyVars_' int2str(iTrials) '.mat'],'bodyMotion','cTimes');
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,length(bodyFiles));
        toc
    end
end
