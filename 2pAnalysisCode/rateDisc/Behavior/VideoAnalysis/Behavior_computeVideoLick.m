function Behavior_computeVideoLick(Animal, Rec, botCam)
% Code to predict lick events in the stimulus peroid based on licks in the decision window.

if ~exist('eyeCam', 'var')
    botCam = 1; %camera that is viewing the face. Default is cam1.
end

Paradigm = 'SpatialDisc';
fPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec '\BehaviorVideo']; %server data path

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end
opts.varCnt = 50; %number of trials used to compute variance map
opts.eyeFrame = 10; %add some pixels to eyeFrame to ensure pupil is properly captured
opts.eyeThresh = .25;
opts.eyeErrode = 2;
opts.maxFrameCnt = 1000; %max number of frames per trial. Earlier frames will be removed from analysis.

% ridge regression variables
ridgeCycles = 5;    %number of cycles for ridge estimation
ridgeSteps = 10;    %number of steps for first cycle
ridgeFolds = 10;     %number of folds for cross validation
maxRidge = 5;       %range of values around ridgeVal that should be tested
ridgeVal = 5;       %initial estimate for ridge parameter
maxFrameCnt = 1000; %max nr of video frames per trial. this should match opts.maxFrameCnt in the Behavior_ConvertSVD code.
dimCnt = 250;       %number of dimensions used for lick analysis


%% load bhv data and check movie files
bhvFile = dir([fileparts(fPath(1:end-1)) filesep Animal '*.mat']);
load([fileparts(fPath(1:end-1)) filesep bhvFile.name]); %load bpod data
movieFiles = dir([fPath '*Video*_' int2str(botCam) '.mj2']);
timeFiles = dir([fPath '*frameTimes*_' int2str(rem(botCam,2)+1) '.mat']); %time stamps for bottom cam

if ~isfield(SessionData,'mouthPos')
    cFile = [fPath movieFiles(1).name];                            
    v = VideoReader(cFile);
    singleFrame = readFrame(v); 
    
    clear v
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(singleFrame); colormap gray; axis image
    mouthPos = ginput(1);
    mouthPos = [round(mouthPos) 50 botCam];
    close(h);
    SessionData.mouthPos = mouthPos;
    save([fileparts(fPath(1:end-1)) filesep bhvFile.name],'SessionData'); %save bpod data
else
    mouthPos = SessionData.mouthPos;
end
mouthInd = [(mouthPos(2)-mouthPos(3)-1) + (1:mouthPos(3)*2+1);(mouthPos(1)-mouthPos(3)-1) + (1:mouthPos(3)*2+1)];

%% loop through trials and perform analysis for pupil, snout and face
Cnt = 0;
imgCnt = zeros(1,3);
leftFrames = cell(1,length(movieFiles));
rightFrames = cell(1,length(movieFiles));

for iTrials = 1:length(movieFiles)
    
    cFile = [fPath movieFiles(iTrials).name];
    v = VideoReader(cFile);
    
    % get timestamps for bottom camera
    load([fPath timeFiles(iTrials).name]);
    if size(frameTimes,1) > opts.maxFrameCnt
        frameTimes = frameTimes(end-opts.maxFrameCnt+1:end);
    end
    botTimes = (frameTimes - SessionData.TrialStartTime(iTrials)) * 86.4*1e3 ;
        
    % find lick frames
    cLickL = {};
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port1In') %check for left licks
        licks = SessionData.RawEvents.Trial{iTrials}.Events.Port1In;
        lickEnd = SessionData.RawEvents.Trial{iTrials}.Events.Port1Out;
        
        licks(licks < SessionData.RawEvents.Trial{iTrials}.States.MoveSpout(end)+0.1) = []; %dont use false licks that occured before spouts were moved in
        licks(licks > max(SessionData.RawEvents.Trial{iTrials}.Events.Port1Out)) = []; %don't use licks that have no following release time
        
        for iLicks = 1:length(licks)
            cLickOn = find(botTimes > licks(iLicks),1); %lick onset time
            cLickEnd = lickEnd(lickEnd - licks(iLicks) > 0); %find correct offset
            cLickEnd = find(botTimes > min(cLickEnd),1); %lick offset time
            cLickL{iLicks} = cLickOn : cLickEnd; %index for current lick
        end
    end
    cLickL = cat(2,cLickL{:});

    cFrames = zeros(v.Height,v.Width,length(cLickL));
    
    
    
    
    rawData = squeeze(importdata([fPath movieFiles(iTrials).name])); %load eye data for current trial    
    rawData = mat2gray(diff(rawData(mouthInd(1,:),mouthInd(2,:),:),1,3));
    
    
    
    rawData = im2uint8(reshape(rawData,size(rawData,1),size(rawData,2),1,[]));
    v = VideoWriter([fPath 'mouthTrace' int2str(iFiles) '.mj2'],'Archival'); %write small mouth video
    open(v);writeVideo(v,rawData);close(v);
    
                                    for iLicks = 1:length(cLickL);
        
        v.CurrentTime = cTimes(cLickL(iLicks));
        cFrames(:,:,iLicks) = readFrame(v);
        
    end
        
    
    
    cLickR = {};
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        licks = SessionData.RawEvents.Trial{iTrials}.Events.Port3In;
        lickEnd = SessionData.RawEvents.Trial{iTrials}.Events.Port3Out;
        licks(licks < SessionData.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        licks(licks > max(SessionData.RawEvents.Trial{iTrials}.Events.Port3Out)) = []; %don't use licks that have no following release time
        
        for iLicks = 1:length(licks)
            cLickOn = find(botTimes > licks(iLicks),1); %lick onset time
            cLickEnd = lickEnd(lickEnd - licks(iLicks) > 0); %find correct offset
            cLickEnd = find(botTimes > min(cLickEnd),1); %lick offset time
            cLickR{iLicks} = cLickOn : cLickEnd; %index for current lick
        end
    end
    cLickR = cat(2,cLickR{:});
    
    
    
    % check for decision period and isolate frames that are non-licking

    if ~isempty(cLickL);
        decTime([cLickL unique(bsxfun(@minus,cLickL(diff([0 cLickL]) > 1) ,(1:3)'))']) = false;  %index for non-licking, also exlude two frames before each lick to avoid confusion
        decTime(unique(bsxfun(@plus,cLickL(diff([cLickL Inf]) > 1) ,(1:3)'))') = false;  %index for non-licking, also exlude two frames after each lick to avoid confusion
    end
    if ~isempty(cLickR);
        decTime([cLickR unique(bsxfun(@minus,cLickR(diff([0 cLickR]) > 1) ,(1:3)'))']) = false;  %index for non-licking, also exlude each frame right before a lick to avoid confusion
        decTime(unique(bsxfun(@minus,cLickR(diff([cLickR Inf]) > 1) ,(1:3)'))') = false;  %index for non-licking, also exlude each frame right before a lick to avoid confusion
    end
    
    % check for waiting period
    waitOn = SessionData.RawEvents.Trial{iTrials}.States.WaitForLever(end);
    waitOff = SessionData.RawEvents.Trial{iTrials}.States.WaitForCam(1);
    waitTime = botTimes > waitOn & botTimes < waitOff; %index for waiting period
    
    % check for stimulus period
    stimOn = SessionData.RawEvents.Trial{iTrials}.States.PlayStimulus(1);
    stimTime = botTimes > stimOn & botTimes < decOn; %index for stimulus period
    
    waitIdx = [waitIdx (find(waitTime)' + botCnt)];
    stimIdx = [stimIdx (find(stimTime)' + botCnt)];
    noLickIdx = [noLickIdx (find(decTime)' + botCnt)];
    lLickIdx = [lLickIdx (cLickL + botCnt)];
    rLickIdx = [rLickIdx (cLickR + botCnt)];    
    clear cLickL cLickR
    
    botCnt = botCnt + length(botTimes);
    
    % save data and increase counter
%     save([fPath 'faceVars_' int2str(iTrials) '.mat'],'eyeVars','snoutMotion','faceMotion','cTimes');
%     Cnt = Cnt + length(cTimes);
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,length(eyeFiles));
        toc
    end
end

%%

minFrames = min([length(noLickIdx) length(waitIdx)]); %minimum framecount in waiting and decision period. 
noLickIdx = randsample(noLickIdx,minFrames); %take an even amount of frames from each period and combine into one index
waitIdx = randsample(waitIdx,minFrames);
noLickIdx = [noLickIdx waitIdx];

%%
% data = [botU(mouthInd,:)*botV(:,lLickIdx) botU(mouthInd,:)*botV(:,noLickIdx) botU(mouthInd,:)*botV(:,rLickIdx)];
data = [botV(:,lLickIdx) botV(:,noLickIdx) botV(:,rLickIdx)];
label = [ones(1,length(lLickIdx)) zeros(1,length(noLickIdx)) -ones(1,length(rLickIdx))];

Mdl = fitcecoc(data',label,'Learners', templateSVM('Standardize',1),...
    'ClassNames',[1 0 -1],'Verbose',2);


stimD = botU* botV(:,stimIdx);
stimD = arrayShrink(stimD,botMask,'split');

lLickD = botU * botV(:,lLickIdx);
lLickD = arrayShrink(lLickD,mouthMask,'split');

rLickD = botU(mouthInd,:) * botV(:,rLickIdx);
rLickD = arrayShrink(rLickD,mouthMask,'split');

decD = botU(mouthInd,:) * botV(:,noLickIdx);
decD = arrayShrink(decD,mouthMask,'split');

        
%% find video components for licking with ridge regression

    trialOn = bhv.TrialStartTime(iTrials) * 86400 + (nanmean([visStim(iTrials) audStim(iTrials)]));
    frameTimes = frameTimes  * 86400 - trialOn;
    trialOn = find(frameTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(frameTimes)));
    trialOn = round(trialOn * (sRate / bhvFrameRate)); %trial on in resampled frames
    
    medianFrameRate(iTrials) = bhvFrameRate;
    peakFrameRate(iTrials) = max(diff(frameTimes*86400));
    

        pupil = mean(eyeVars.axes,2); %pupil diameter
        idx = zscore(pupil) < -2; %index for blinks
        pupil(idx) = mean(pupil); %replace blinks with pupil average
        pupil = [repmat(pupil(1),21,1); pupil; repmat(pupil(end),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        pupil = resample(pupil, sRate, bhvFrameRate); %resample to match imaging data framerate
        offset = ceil(21 * sRate / bhvFrameRate); %required offset to remove padding after resampling
        pupil = pupil(offset +1 : end - offset); %remove padds
        pupil = smooth(pupil,'rlowess'); %do some smoothing
        pupilR{iTrials} = single(pupil(trialOn:trialOn + (size(timeR,1) - 1))); %only use trial-relevant frames
    
        



if exist([cPath 'lickRidgeTest.mat']) == 2
    load([cPath 'lickRidgeTest.mat'])
else    
    ridgeRange = (-maxRidge : maxRidge / (ridgeSteps-1) * 2 : maxRidge) /2 + ridgeVal; %range of values to be tested
    ridgeRange(ridgeRange <= 0) = [];
    
    if ~isunix; figure; hold on; end
    for iCycles = 1:ridgeCycles
        ridgeRMSE = NaN(ridgeFolds,length(ridgeRange));
        Cnt = 0;
        for iSteps = ridgeRange
            Cnt = Cnt + 1;
            foldCnt = floor(size(V,2) / ridgeFolds);
            
            for iFolds = 1:ridgeFolds
                
                cIdx = true(1,size(V,2));
                cIdx(((iFolds - 1)*foldCnt) + (1:foldCnt)) = false;
                
                ridgeR = [fullR(cIdx,:); diag(ones(1,size(fullR,2))) * iSteps];
                ridgeBeta = (ridgeR' * ridgeR) \ ridgeR' * [V(:,cIdx) zeros(size(V,1),size(ridgeR,2))]'; %compute beta weight maps with ridge penalty for selected subset
                temp = (V(:,~cIdx)' - (fullR(~cIdx,:) * ridgeBeta)).^2; %compute error when predicting the remaining dataset
                
                ridgeRMSE(iFolds,Cnt) = gather(sqrt(mean(temp(:)))); %compute error in remaining dataset
                fprintf('iFold: %d/%d, iSteps: %d/%d, iCycles: %d/%d\n', iFolds, ridgeFolds, Cnt, length(ridgeRange), iCycles, ridgeCycles);
                
                if iFolds == ridgeFolds && iSteps == ridgeRange(end)
                    
                    ridgeTest{iCycles}{1} = ridgeRange;
                    ridgeTest{iCycles}{2} = ridgeRMSE;
                    if ~isunix; plot(ridgeRange, nanmean(ridgeRMSE),'o-','linewidth',2); axis square; end
                    
                    if iCycles ~= ridgeCycles
                        
                        [~,temp] = min(nanmean(ridgeRMSE)); %find ridge value that gave lowest error
                        if temp ~= size(ridgeRMSE,2)
                            maxRidge = ceil(maxRidge / 2); %decrease the range for next run
                        end
                        ridgeRange = (-maxRidge : maxRidge / (ridgeSteps-1) * 2 : maxRidge) /2 + ridgeRange(temp); %range of values to be tested
                        ridgeRange(ridgeRange <= 0) = [];
                        
                    end
                end
            end
        end
    end
    
    [~,temp] = min(nanmean(ridgeRMSE)); %find ridge value that gave lowest error
    ridgeVal = ridgeRange(temp); %use this value for regression on full dataset
    save([cPath 'lickRidgeTest.mat'], 'ridgeTest', 'ridgeVal')
end



