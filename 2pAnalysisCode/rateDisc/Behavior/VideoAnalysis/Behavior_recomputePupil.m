function Behavior_recomputePupil(Animal, Rec, eyeCam, newEye)
% Code to convert face variables from compressed movies. Will compute face
% motion, body motion and pupil diameter and save in faceVars.mat.
% This code is meant to be used after converting videos to SVD. 
% Needs eyeTrace.mj2 videos for pupil and V and U components for face variables.
% If newEye = true, will compute eyeTrace from SVD data

if ~exist('eyeCam', 'var')
    eyeCam = 1; %camera that is viewing the face. Default is cam1.
end

if ~exist('newEye', 'var')
    newEye = false; %camera that is viewing the face. Default is cam1.
end

Paradigm = 'SpatialDisc';
fPath = ['Y:\data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec '\BehaviorVideo']; %server data path

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end


%% check files and load SVD data
load([fPath '\bhvOpts.mat'])
load([fPath '\segInd'  num2str(eyeCam) '.mat'])

Frame = Behavior_rebuildSegmentFrame(fPath, eyeCam);





load([fPath '\rawMean_Cam' num2str(eyeCam) '.mat'])
load([fPath '\SVD_Cam ' num2str(eyeCam) '.mat'])

bhvFile = dir([fileparts(fPath(1:end-1)) '\' Animal '*.mat']);
load([fileparts(fPath(1:end-1)) '\' bhvFile.name]); %load bpod data
snoutPos = SessionData.snoutPos / 2; %get snout position and camera for face
eyeFiles = dir([fPath '*eyeTrace*.mj2']);

U = U((snoutPos(2)-snoutPos(3)-1) + (1:snoutPos(3)*2+1), (snoutPos(1)-snoutPos(3)-1) + (1:snoutPos(3)*2+1),:); %spatial components for snout motion
U = reshape(U,size(U,1)*size(U,2),[]);
        
%% loop through trials and perform analysis
tic;
Cnt = 0;
for iTrials = 1:length(eyeFiles)
    
    load([fPath 'faceVars_' int2str(iTrials) '.mat']); %load last dataset. This is mainly to get cTimes which can be saved down again later.

    eyeTrace = squeeze(importdata([fPath 'eyeTrace_' int2str(iTrials) '.mj2'])); %load eye data for current trial
    eyeTrace = mat2gray(eyeTrace);
    
    eyeVars.orient = zeros(size(eyeTrace,3),1);
    eyeVars.center = zeros(size(eyeTrace,3),2);
    eyeVars.axes = zeros(size(eyeTrace,3),2);
    eyeVars.solidity = zeros(size(eyeTrace,3),1);
    
    for iFrames = 1:size(eyeTrace,3)
        temp = Behavior_EyeCheck(eyeTrace(:,:,iFrames),0.25,0); %get pupil data
        eyeVars.center(iFrames,:) = temp.center; %circle center
        eyeVars.axes(iFrames,:) = temp.axes; %circle axis
        eyeVars.orient(iFrames) = temp.orient; %circle orientation
        eyeVars.solidity(iFrames) = temp.solidity; %circle solidity
    end
    
    % compute snout motion
    snoutMotion = U * V(:,Cnt + (1:size(eyeTrace,3)));
    snoutMotion = abs(diff(snoutMotion,1,2));
    snoutMotion = mean(snoutMotion);
    snoutMotion = [snoutMotion(1) snoutMotion];
    
    %coompute face motion
    faceMotion = abs(diff(V(:,Cnt + (1:size(eyeTrace,3))),1,2));
    faceMotion = mean(faceMotion);
    faceMotion = [faceMotion(1) faceMotion];

    % save data and increase counter
    save([fPath 'faceVars_' int2str(iTrials) '.mat'],'eyeVars','snoutMotion','faceMotion','cTimes');
    Cnt = Cnt + size(eyeTrace,3);
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,length(eyeFiles));
        toc
    end
end