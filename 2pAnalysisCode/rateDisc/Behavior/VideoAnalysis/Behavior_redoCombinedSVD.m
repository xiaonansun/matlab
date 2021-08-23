function Behavior_redoCombinedSVD(cPath)
% code to redo final SVD for behavioral videos segments
tic;
load([cPath 'bhvOpts.mat']); %get settings for previous SVD
bhvOpts.trialDur = bhvOpts.framesPerTrial / bhvOpts.targRate; %amount of frames per trial in trial-aligned video data
nrCams = 2;

fPath = [fileparts(cPath(1:end-1)) filesep];
bhvFile = strsplit(fPath,filesep);
bhvFile = dir([fPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
load([fPath bhvFile.name]); %load behavior data

%% load svd data and merge different segments together
allV = cell(1,nrCams);
allTimes = cell(1,nrCams);
for iCams = 1:nrCams %get segments and frameTimes for each cam
    load([cPath 'segInd' int2str(iCams) '.mat'],'ind'); %get index for reconstruction
    load([cPath 'SVD_Cam' int2str(iCams) '-frameTimes.mat'],'totalFrameTimes','frameCnt');
    
    % check if timestamps are shifted by an hour. Apparently that can happen.
    timeCheck = (SessionData.TrialStartTime(1)*86400) - (totalFrameTimes(1) * 86400); %time difference between first acquired frame and onset of first trial
    if timeCheck > 3590 && timeCheck < 3610 %timeshift by one hour (+- 10seconds)
        warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
        totalFrameTimes = totalFrameTimes + (1/24); %add one hour
    elseif timeCheck > 30
        error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
    end
    
    allTimes{iCams} = totalFrameTimes; %get timestamps for current cam
    allFrameCnt{iCams} = frameCnt;
    bhvOpts.frameRate(iCams) = round(1/median(diff(allTimes{iCams} * 86400))); %determine frameRate for current cam
    
    % load compressed video segments
    allV{iCams} = cell(1,size(ind,2));
    for iSegs = 1 : size(ind,2)
        data = load([cPath 'SVD_Cam' int2str(iCams) '-Seg' int2str(iSegs) '.mat'],'V');
        allV{iCams}{iSegs} = data.V; clear data
    end
    allV{iCams} = cat(2,allV{iCams}{:}); % combine in on larger array
end

%% go through each trial and build combined dataset for each cam
for x = 1 : nrCams
    alignV{x} = cell(1,size(SessionData.RawEvents.Trial,2));
end

for iTrials = 1:size(SessionData.RawEvents.Trial,2)
    try
        stimTime(iTrials) = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    catch
        stimTime(iTrials) = NaN;
    end
    
    for iCams = 1:nrCams
        try
            trialOn = (SessionData.TrialStartTime(iTrials) * 86400) + (stimTime(iTrials) - bhvOpts.preStimDur); % trial onset times
            cTimes = allTimes{iCams} * 86400 - trialOn; % get frame times for current cam
            trialOn = find(cTimes > 0,1); %trial on in frames
            trialOff = find(cTimes > bhvOpts.trialDur,1); %trial off in frames
            
            trialTime = cTimes(trialOn:trialOff) - cTimes(trialOn);
            idx = trialTime < bhvOpts.trialDur; %don't use late frames
            trialTime = trialTime(idx);
            timeLeft = bhvOpts.trialDur - trialTime(end); %check if there is missing time at the end of a trial
            
            cData = allV{iCams}(trialOn:trialOff,:); %get data based on time vector matching duration of a trial
            cData = cData(idx, :); %don't use late frames
            
            if (timeLeft < bhvOpts.trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
                addTime = trialTime(end) + (1/bhvOpts.targRate : 1/bhvOpts.targRate : timeLeft + 1/bhvOpts.targRate); %add some dummy times to make complete trial
                trialTime = [trialTime' addTime];
                cData = [cData; repmat(cData(1,:),length(addTime),1)];
            end
            
            if any(diff(trialTime) > 0.5) %dont use current trial if frames are no more than 0.5s apart
                error;
            end
            
            %create time vector to capture potential irregular sampling
            timePad = 1/bhvOpts.targRate:1/bhvOpts.targRate:1/bhvOpts.targRate*21;
            trialTime = trialTime'+timePad(end)+1/bhvOpts.frameRate(iCams);
            trialTime = [timePad trialTime' timePad+trialTime(end)];
            cData = [repmat(cData(1,:),21,1); cData; repmat(cData(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
            cData = resample(double(cData), trialTime, bhvOpts.targRate); %resample to set target framerate
            cData = conv2(cData',ones(1,bhvOpts.smth)/bhvOpts.smth,'same')'; %smooth trace with moving average of 'smth' points
            alignV{iCams}{iTrials} = cData(22 : 21 + bhvOpts.framesPerTrial,:); %remove padds and select trial-relevant frames
        catch
            alignV{iCams}{iTrials} = NaN(bhvOpts.framesPerTrial,size(allV{iCams},2),'single'); %dont use this trial
        end
    end
    if rem(iTrials,50) == 0
        fprintf(1, 'Current run is %d out of %d\n', iTrials, size(SessionData.RawEvents.Trial,2));
        toc
    end
end

% combined V segments from all cams
for x = 1 : nrCams
    alignV{x} = cat(1,alignV{x}{:}); %combined V matrix
end
alignV = cat(2,alignV{:}); %combined V matrix
camIdx = ~isnan(mean(alignV,2)); %index for trials where bhv frames were available

% create second SVD
tic; disp('Compute combined SVD.');
[V,S,vidU]  = svd(alignV(camIdx,:), 'econ'); clear alignV allV
V = V * S;
V = V(:,1:bhvOpts.nSVD); %get request nr of dimensions
vidU = vidU(:,1:bhvOpts.nSVD)'; %the resulting U is nSVD x regs

%zscore video regressor, save mean and std for later reconstruction
meanV = mean(V);
stdV = std(V);
V = bsxfun(@minus,V,meanV);
V = bsxfun(@rdivide,V,stdV);

% rebuild
vidV = NaN(length(camIdx),bhvOpts.nSVD,'single');
vidV(camIdx,:) = V; clear V
save([cPath 'SVD_CombinedSegments.mat'], 'vidV', 'vidU', 'meanV' ,'stdV','allTimes');
save([cPath 'bhvOpts.mat'], 'bhvOpts');
disp('Done.'); toc;