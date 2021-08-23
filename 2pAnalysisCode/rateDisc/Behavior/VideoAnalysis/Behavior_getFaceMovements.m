function [eyeVars, snoutMotion, noseMotion] = Behavior_getFaceMovements(cVideo, fPath, eyeThresh, repairOutliers)
% code to extract facial features from a single video file. will isolate
% and save eyeTrace and compute pupil diameter + snout and nose motion.

if ~exist('eyeThresh', 'var') || isempty(eyeThresh) % threshold when computing pupil diameter
    eyeThresh(1) = 0;
    eyeThresh(2) = 0.075;
end

if ~exist('repairOutliers', 'var') % flag to repair outliers in pupil estimates. This can take quite a long time - make sure rejection parameters are appropriate for your data.
    repairOutliers = false;
end
frameChange = 1; %maximum change in diameter estimate from one frame to the next (pixels)
slowDiff = 5; %maximum difference between estimate and a smoothed estimate (pixels)
    
% check for behavioral file containing position for eye, snout and nose
if ~exist([fPath 'SessionData.mat'], 'file')
    SessionData = Behavior_addNose(fPath,[],1); %add eye, snout and nose positions
    SessionData.eyePos(1:3) = round(SessionData.eyePos(1:3) / 2);
    SessionData.snoutPos(1:3) = round(SessionData.snoutPos(1:3) / 2);
    SessionData.nosePos(1:3) = round(SessionData.nosePos(1:3) / 2);
    save([fPath 'SessionData.mat'], 'SessionData');
else
    load([fPath 'SessionData.mat']);
end
eyePos = SessionData.eyePos;
snoutPos = SessionData.snoutPos;
nosePos = SessionData.nosePos;

% get eye movie and compute pupil diameter
eyeTrace = cVideo((eyePos(2)-eyePos(3)-1) + (1:eyePos(3)*2+1), (eyePos(1)-eyePos(3)-1) + (1:eyePos(3)*2+1), :);
eyeTrace = mat2gray(eyeTrace);
eyeTrace = im2uint8(reshape(eyeTrace,size(eyeTrace,1),size(eyeTrace,2),1,[]));
eyeTrace = padarray(eyeTrace, 8 - [rem(size(eyeTrace,1),8), rem(size(eyeTrace,2),8)], 0, 'post'); %pad array to be dividable by 8
eyeTrace = reshape(eyeTrace, size(eyeTrace,1),size(eyeTrace,2),1,[]);
        
v = VideoWriter([fPath 'eyeTrace.mp4'], 'MPEG-4'); %write small eye video
open(v); writeVideo(v,eyeTrace); close(v);

% get pupil diameter
eyeTrace = squeeze(mat2gray(eyeTrace));
eyeVars.orient = zeros(size(eyeTrace,3),1);
eyeVars.center = zeros(size(eyeTrace,3),2);
eyeVars.axes = zeros(size(eyeTrace,3),2);
eyeVars.solidity = zeros(size(eyeTrace,3),1);

for iFrames = 1:size(eyeTrace,3)
    temp = Behavior_EyeCheck(eyeTrace(:,:,iFrames),eyeThresh(2),eyeThresh(1),false); %get pupil data
    eyeVars.center(iFrames,:) = temp.center; %circle center
    eyeVars.axes(iFrames,:) = temp.axes; %circle axis
    eyeVars.orient(iFrames) = temp.orient; %circle orientation
    eyeVars.solidity(iFrames) = temp.solidity; %circle solidity
    if rem(iFrames,1000) == 0
        fprintf('Pupil estimate: %g/%g frames complete.\n', iFrames,size(eyeTrace,3));
    end
end

% find bad pupil estimates and try to repair.
if repairOutliers
    trace = mean(eyeVars.axes,2);
    if sum(trace == 0) ~= length(trace)
        trace(trace == 0) = NaN; %replace zero-estimates with NaNs
        trace = fillgaps(trace, 100); %fill NaN throgh interpolation
        smoothTrace = smooth(trace,10,'rlowess'); %use smooth pupil trace as reference for outliers
    else
        smoothTrace = trace;
    end
    smoothDiff = abs(trace-smoothTrace); %difference from smoothed pupil trace
    smoothStd = std(smoothDiff);
    smoothDiff = smoothDiff ./ smoothStd; %put to STDs
    idx = find([0;diff(zscore(trace))] < -frameChange | [0;diff(zscore(trace))] > frameChange | smoothDiff > slowDiff); %potentially bad results. Try to do better by finding result that is closest to average.
    
    for iFrames = 1:length(idx)
        clear temp tDiff
        Cnt = 0;
        for x = [-6 -2 2 6] %try eroding or fusing patches to get better pupil estimate
            Cnt = Cnt + 1;
            temp{Cnt} = Behavior_EyeCheck(eyeTrace(:,:,idx(iFrames)),eyeThresh(2),eyeThresh(1),false,true,x); %get pupil data. Force individual areas to fuse together.
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
end
        
% compute snout motion
idx1 = (snoutPos(2)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
idx2 = (snoutPos(1)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
if any(idx1 > size(cVideo,1)) || any(idx2 > size(cVideo,2))
    if iFiles == 1
        warning('Snout index is too large. Clipped to fit into image.');
    end
    idx1(idx1 > size(cVideo,1)) = [];
    idx2(idx2 > size(cVideo,2)) = [];
end
snoutMotion = cVideo(idx1, idx2, :);
snoutMotion = diff(snoutMotion,1,3);
snoutMotion = reshape(snoutMotion,[],size(cVideo,3)-1);
snoutMotion = mean(snoutMotion);
snoutMotion = [snoutMotion(1) snoutMotion]';

%compute nose motion
idx1 = (nosePos(2)-nosePos(3)-1) + (1:nosePos(3)*2+1);
idx2 = (nosePos(1)-nosePos(3)-1) + (1:nosePos(3)*2+1);
if any(idx1 > size(cVideo,1)) || any(idx2 > size(cVideo,2))
    if iFiles == 1
        warning('Nose index is too large. Clipped to fit into image.');
    end
    idx1(idx1 > size(cVideo,1)) = [];
    idx2(idx2 > size(cVideo,2)) = [];
end
noseMotion = cVideo(idx1, idx2, :);
noseMotion = diff(noseMotion,1,3);
noseMotion = reshape(noseMotion,[],size(cVideo,3)-1);
noseMotion = mean(noseMotion);
noseMotion = [noseMotion(1) noseMotion]';
