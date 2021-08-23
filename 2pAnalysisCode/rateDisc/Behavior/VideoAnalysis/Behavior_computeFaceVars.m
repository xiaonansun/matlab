function Behavior_computeFaceVars(cPath, Animal, Rec, eyeCam, reAssign)
% Code to convert face variables from compressed movies. Will compute face
% motion, body motion and pupil diameter and save in faceVars.mat.
% This code is meant to be used after converting videos to SVD. 
% Needs eyeTrace.mj2 videos for pupil and V and U components for face variables.

if ~exist('eyeCam', 'var')
    eyeCam = 1; %camera that is viewing the face. Default is cam1.
end

if ~exist('reAssign', 'var')
    reAssign = false; %flag that snout and nose position should be re-assigned
end

Paradigm = 'SpatialDisc';
fPath = [cPath Animal filesep Paradigm filesep Rec '\BehaviorVideo']; %server data path

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end
noseSize = 20; %half-size of nose rectangle if reassigned

if exist([fPath 'bhvOpts.mat'],'file') %check bhvOpts and load if present
    load([fPath 'bhvOpts']);
else
    bhvOpts.maxFrameCnt = 500; %max number of frames per trial. Earlier frames will be removed from analysis.
end

%% check files and load raw data
eyeFiles = dir([fPath '*eyeTrace*.mj2']);
faceMovies = dir([fPath Animal '*_' num2str(eyeCam) '.mj2']);
bodyMovies = dir([fPath Animal '*_' num2str(rem(eyeCam,2)+1) '.mj2']);
bhvFile = dir([fileparts(fPath(1:end-1)) '\' Animal '_' Paradigm '*.mat']);
load([fileparts(fPath(1:end-1)) '\' bhvFile.name]); %load bpod data

%% check if nose position is assigned or nose/snout position should be re-assigned
if ~isfield(SessionData, 'nosePos') || reAssign
    SessionData.nosePos(3) = noseSize; 
    SessionData.nosePos(4) = eyeCam; 
    h1 = figure;
    v = VideoReader([fPath faceMovies(1).name]);
    temp = readFrame(v);imgSize = size(temp);clear v
    imagesc(temp);axis image; colormap gray
    text(round(imgSize(1) * 0.25), imgSize(2) - round(imgSize(2)/1.1),'click SNOUT','FontSize',30,'color','b')
    [xPos,yPos] = ginput(1);
    SessionData.snoutPos(1) = round(xPos);
    SessionData.snoutPos(2) = round(yPos);
    
    h2 = figure;
    ax = subplot(1,2,1);
    imagesc(temp);axis image; colormap gray
    Frame = [SessionData.snoutPos(1)-SessionData.snoutPos(3) SessionData.snoutPos(2)-SessionData.snoutPos(3) SessionData.snoutPos(3)*2 SessionData.snoutPos(3)*2];
    rectangle('Position',Frame,'linewidth',2,'edgecolor','b','parent',ax);
    title('Snout position')

    figure(h1);
    v = VideoReader([fPath faceMovies(1).name]);
    temp = readFrame(v);imgSize = size(temp); clear v
    imagesc(temp);axis image; colormap gray
    text(round(imgSize(1) * 0.25), imgSize(2) - round(imgSize(2)/1.1),'click NOSE','FontSize',30,'color','r')
    [xPos,yPos] = ginput(1);
    SessionData.nosePos(1) = round(xPos);
    SessionData.nosePos(2) = round(yPos);
    
    figure(h2);
    ax = subplot(1,2,2);
    imagesc(temp);axis image; colormap gray
    Frame = [SessionData.nosePos(1)-SessionData.nosePos(3) SessionData.nosePos(2)-SessionData.nosePos(3) SessionData.nosePos(3)*2 SessionData.nosePos(3)*2];
    rectangle('Position',Frame,'linewidth',2,'edgecolor','r','parent',ax);
    title('Nose position');
    
    %save bpod data
    save([fileparts(fPath(1:end-1)) '\' bhvFile.name], 'SessionData');
    close(h1); drawnow;
end
snoutPos = SessionData.snoutPos; %get snout position and camera for face
nosePos = SessionData.nosePos; %get snout position and camera for face

% loop through trials and perform analysis
tic;
for iTrials = 1:length(eyeFiles)
    
    cFile = [fPath faceMovies(iTrials).name];
    rawData = squeeze(importdata(cFile));
    if size(rawData,3) > bhvOpts.maxFrameCnt
        rawData = rawData(:,:,end-bhvOpts.maxFrameCnt+1:end);
    end
    
    % compute snout motion
    idx1 = (snoutPos(2)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
    idx2 = (snoutPos(1)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
    if any(idx1 > size(rawData,1)) || any(idx2 > size(rawData,2))
        if iTrials == 1
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
    
    %compute face motion
    faceMotion = diff(rawData,1,3);
    faceMotion = reshape(faceMotion,[],size(rawData,3)-1);
    faceMotion = mean(faceMotion);
    faceMotion = [faceMotion(1) faceMotion];
        
    %compute nose motion
    idx1 = (nosePos(2)-nosePos(3)-1) + (1:nosePos(3)*2+1);
    idx2 = (nosePos(1)-nosePos(3)-1) + (1:nosePos(3)*2+1);
    if any(idx1 > size(rawData,1)) || any(idx2 > size(rawData,2))
        if iTrials == 1
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

    % compute body motion
    cFile = [fPath bodyMovies(iTrials).name];
    rawData = squeeze(importdata(cFile));
    if size(rawData,3) > bhvOpts.maxFrameCnt
        rawData = rawData(:,:,end-bhvOpts.maxFrameCnt+1:end);
    end
    bodyMotion = diff(rawData,1,3);
    bodyMotion = reshape(bodyMotion,[],size(rawData,3)-1);
    bodyMotion = mean(bodyMotion);
    bodyMotion = [bodyMotion(1) bodyMotion];
    
    % save data and increase counter
    save([fPath 'faceVars_' int2str(iTrials) '.mat'],'snoutMotion','noseMotion','faceMotion','bodyMotion','-append');
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,length(eyeFiles));
        toc
    end
end
