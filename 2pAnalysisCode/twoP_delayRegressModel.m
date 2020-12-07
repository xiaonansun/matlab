function twoP_delayRegressModel(cPath,Animal,Rec)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%% general variables
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sRate = round(30.8987); % Sampling rate of imaging in Hz
preStimDur = ceil(1.8*sRate) / sRate; % Duration of trial before lever grab in seconds
postStimDur = ceil(4.5*sRate) / sRate; % Duration of trial after lever grab onset in seconds
frames = round((preStimDur + postStimDur) * sRate); %nr of frames per trial
trialDur = (frames * (1/sRate)); %duration of trial in seconds

%other variables
mPreTime = ceil(0.5*sRate) / sRate;  % precede motor events to capture preparatory activity in seconds
mPostTime = ceil(1*sRate) / sRate;   % follow motor events for mPostStim in seconds
motorIdx = [-((mPreTime * sRate): -1 : 1) 0 (1:(mPostTime * sRate))]; %index for design matrix to cover pre- and post motor action
tapDur = 0.25;      % minimum time of lever contact, required to count as a proper grab.
noWhisk = 5;        % minimum nrs of frames before and after a lick to check for whisking. Whisks should not be counted during licking to reduce confusion.
piezoLine = 5;      % channel in the analog data that contains data from piezo sensor
aCenter = 3001;     % avg. stimulus onset in analog data

bhvDimCnt = 50;    % number of dimensions from behavioral videos that are used as regressors.
gaussShift = 1;     % inter-frame interval between regressors. Will use only every 'gaussShift' regressor and convolve with gaussian of according FHWM to reduce total number of used regressors.

trialSegments = {1:54 55:81 89:132 140:162 170:188}; %trial segments to determine motion.
% motor labels are used to orthogonalize the motor regressors in the model against the trial-time

%% load data
load([cPath 'data']); %load 2p data
bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data
SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds

% ensure there are not too many trials in the dataset
bTrials = data.trialNumbers;
bTrials(~ismember(data.trialNumbers,data.bhvTrials)) = []; %don't use trials that have problems with trial onset times
bTrials(SessionData.DidNotChoose(bTrials) | SessionData.DidNotLever(bTrials) | ~SessionData.Assisted(bTrials)) = []; %don't use unperformed/assisted trials

data.dFOF(:,:,~ismember(data.trialNumbers,bTrials)) = [];
data.analog(:,:,~ismember(data.trialNumbers,bTrials)) = [];
bhv = selectBehaviorTrials(SessionData,bTrials); %only use selected trials
trialCnt = length(bTrials);

%% load behavior data
if exist([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'file') ~= 2 %check if svd behavior exists on hdd and pull from server otherwise
    if exist([cPath 'BehaviorVideo' filesep]) ~= 2
        mkdir([cPath 'BehaviorVideo' filesep]);
    end
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'FilteredPupil.mat'],[cPath 'BehaviorVideo' filesep 'FilteredPupil.mat']);
end

load([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'vidV'); %load behavior video data
V1 = vidV(:,1:bhvDimCnt); clear vidV %behavioral video regressors

load([cPath 'BehaviorVideo' filesep 'FilteredPupil.mat']); %load pupil data
%check if timestamps from pupil data are shifted against bhv data
timeCheck = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
if timeCheck > 3590 && timeCheck < 3610 %timeshift by one hour (+- 10seconds)
    warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
    for iTrials = 1 : length(pTime)
        pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
    end
elseif timeCheck > 30 || timeCheck < -30
    error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
end
    
%% find events in BPod time - All timestamps are relative to stimulus onset event to synchronize to imaging data later
% pre-allocate vectors
lickL = cell(1,trialCnt);
lickR = cell(1,trialCnt);
leverIn = NaN(1,trialCnt);
levGrabL = cell(1,trialCnt);
levGrabR = cell(1,trialCnt);
levReleaseL = cell(1,trialCnt);
levReleaseR = cell(1,trialCnt);
water = NaN(1,trialCnt);

for iTrials = 1:trialCnt
        
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    
    try
        stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab(iTrials); %time of stimulus onset - measured from soundcard
    catch
        stimTime(iTrials) = NaN;
    end
     
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lickL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port1In;
        lickL{iTrials}(lickL{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickL{iTrials} = lickL{iTrials} - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        lickR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port3In;
        lickR{iTrials}(lickR{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickR{iTrials} = lickR{iTrials} - stimGrab(iTrials);
    end
    
    leverIn(iTrials) = min(bhv.RawEvents.Trial{iTrials}.States.Reset(:)) - stimGrab(iTrials); %first reset state causes lever to move in
        
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimGrab(iTrials);
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimGrab(iTrials);
    end
    
    if ~isnan(bhv.RawEvents.Trial{iTrials}.States.Reward(1)) %check for reward state
        water(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
    end
end

%normalize response side between -1 and 1
bhv.ResponseSide(bhv.ResponseSide == 1) = -1;
bhv.ResponseSide(bhv.ResponseSide == 2) = 1;
SessionData.ResponseSide(SessionData.ResponseSide == 1) = -1;
SessionData.ResponseSide(SessionData.ResponseSide == 2) = 1;
maxStimRegs = length(min(round((preStimDur + stimTime) * sRate)) : (preStimDur + postStimDur) * sRate); %maximal number of required stimulus regressors

%% build regressors - create design matrix based on event times
%basic time regressors
timeR = logical(diag(ones(1,frames)));

lGrabR = cell(1,trialCnt);
lGrabRelR = cell(1,trialCnt);
rGrabR = cell(1,trialCnt);
rGrabRelR = cell(1,trialCnt);
lLickR = cell(1,trialCnt);
rLickR = cell(1,trialCnt);
leverInR = cell(1,trialCnt);

lVisStimR = cell(1,trialCnt);
rVisStimR = cell(1,trialCnt);
lAudStimR = cell(1,trialCnt);
rAudStimR = cell(1,trialCnt);

visRewardR = cell(1,trialCnt);
audRewardR = cell(1,trialCnt);
prevRewardR = cell(1,trialCnt);

visChoiceR = cell(1,trialCnt);
audChoiceR = cell(1,trialCnt);
prevChoiceR = cell(1,trialCnt);
prevModR = cell(1,trialCnt);

waterR = cell(1,trialCnt);
fastPupilR = cell(1,trialCnt);
slowPupilR = cell(1,trialCnt);

whiskR = cell(1,trialCnt);
noseR = cell(1,trialCnt);
piezoR = cell(1,trialCnt);
%%
tic
for iTrials = 1:trialCnt
    %% vis/aud stim - regressors cover the remaining trial after stimulus onset
    stimIdx = round((preStimDur + stimTime(iTrials)) * sRate) : round((preStimDur + postStimDur)*sRate); %index for which part of the trial should be covered by stim regressors
       
    % vision
    lVisStimR{iTrials} = false(frames, maxStimRegs);
    rVisStimR{iTrials} = false(frames, maxStimRegs);
    if bhv.StimType(iTrials) == 1 || bhv.StimType(iTrials) == 3 %visual or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lVisStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        else
            rVisStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        end
    end
    
    % audio
    lAudStimR{iTrials} = false(frames, maxStimRegs);
    rAudStimR{iTrials} = false(frames, maxStimRegs);
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 3 %auditory or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lAudStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        else
            rAudStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        lVisStimR{iTrials} = lVisStimR{iTrials}(:,1:gaussShift:end);
        rVisStimR{iTrials} = rVisStimR{iTrials}(:,1:gaussShift:end);
        lAudStimR{iTrials} = lAudStimR{iTrials}(:,1:gaussShift:end);
        rAudStimR{iTrials} = rAudStimR{iTrials}(:,1:gaussShift:end);
    end
        
    %% lick regressors
    lLickR{iTrials} = false(frames, length(motorIdx));
    rLickR{iTrials} = false(frames, length(motorIdx));
    
    for iRegs = 0 : length(motorIdx)-1
        licks = lickL{iTrials} - (mPreTime - (iRegs * 1/sRate));
        lLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        licks = lickR{iTrials} - (mPreTime - (iRegs * 1/sRate));
        rLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    if gaussShift > 1
        % subsample regressors
        lLickR{iTrials} = lLickR{iTrials}(:,1:gaussShift:end);
        rLickR{iTrials} = rLickR{iTrials}(:,1:gaussShift:end);
    end      
    
    %% lever in
    leverInR{iTrials} = false(frames,frames);
%     leverShift = round((preStimDur + leverIn(iTrials))* sRate); %timepoint in frames when lever moved in, relative to lever grab
%     
%     if ~isnan(leverShift)
%         if leverShift > 0 %lever moved in during the recorded trial
%             leverInR{iTrials}(:, 1: frames - leverShift) = timeR(:, leverShift+1:end);
%         else %lever was present before data was recorded
%             leverInR{iTrials}(:, abs(leverShift) + 1 : end) = timeR(:, 1: frames + leverShift);
%         end
%     end
    
    if gaussShift > 1
        % subsample regressors
        leverInR{iTrials} = leverInR{iTrials}(:,1:gaussShift:end);
    end

    %% choice and reward
    stimShift = round((stimTime(iTrials))* sRate - sRate); %timepoint in frames when the stimulus was presented. This is the shift relative to the expection that the stimulus comes up 1s after grabing the lever.
    
    stimTemp = false(frames,frames);
    if stimShift > 0 %stim came later than 1s from lever grab
        stimTemp(:, 1: frames - stimShift) = timeR(:, stimShift+1:end);
    else %stim came earlier than 1s from lever grab
        stimTemp(:, abs(stimShift) + 1 : end) = timeR(:, 1: frames + stimShift);
    end
    stimTemp(:,end-4:end) = []; %don't use the last 5 timepoints to avoid rank defficient design matrix

    visRewardR{iTrials} = single((stimTemp .* (bhv.StimType(iTrials)  == 1)) * (-1 + (2 * bhv.Rewarded(iTrials)))); %visual trial, 1 for reward, -1 for no reward
    audRewardR{iTrials} = single((stimTemp .* (bhv.StimType(iTrials)  == 2)) * (-1 + (2 * bhv.Rewarded(iTrials)))); %audio trial, 1 for reward, -1 for no reward
    
    visChoiceR{iTrials} = single(stimTemp .* bhv.ResponseSide(iTrials) .* (bhv.StimType(iTrials)  == 1)); %visual choice regressor, takes -1 for left and 1 for right responses
    audChoiceR{iTrials} = single(stimTemp .* bhv.ResponseSide(iTrials) .* (bhv.StimType(iTrials)  == 2)); %auditory choice regressor, takes -1 for left and 1 for right responses
       
    if iTrials == 1 %don't use first trial
        prevRewardR{iTrials} = NaN(size(stimTemp));
        prevChoiceR{iTrials} = NaN(size(stimTemp));
        prevModR{iTrials} = NaN(size(stimTemp));
        
    else %for all subsequent trials
        prevChoiceR{iTrials} = single(stimTemp .* SessionData.ResponseSide(bTrials(iTrials)-1)); % check if previous trial was rewarded (1), punished (-1). Goes to 0 if not performed        
        
        if SessionData.StimType(bTrials(iTrials)-1) == 1 || SessionData.StimType(bTrials(iTrials)-1) == 3
            prevModR{iTrials} = single(stimTemp); % if previous trial was vision (1)
        elseif SessionData.StimType(bTrials(iTrials)-1) == 2 || SessionData.StimType(bTrials(iTrials)-1) == 3
            prevModR{iTrials} = single(stimTemp) * -1; % if previous trial was audio (-1)
        end
            
        if ~SessionData.Rewarded(bTrials(iTrials)-1) %last trial was not rewarded
            if SessionData.DidNotChoose(bTrials(iTrials)-1) %no response in last trial
                prevChoiceR{iTrials} = zeros(size(stimTemp),'single');
                prevRewardR{iTrials} = zeros(size(stimTemp),'single');
            else % punished
                prevRewardR{iTrials} = single(stimTemp .* -1);
            end
        else
            prevRewardR{iTrials} = single(stimTemp);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        visRewardR{iTrials} = visRewardR{iTrials}(:,1:gaussShift:end);
        audRewardR{iTrials} = audRewardR{iTrials}(:,1:gaussShift:end);
        prevRewardR{iTrials} = prevRewardR{iTrials}(:,1:gaussShift:end);
        visChoiceR{iTrials} = visChoiceR{iTrials}(:,1:gaussShift:end);
        audChoiceR{iTrials} = audChoiceR{iTrials}(:,1:gaussShift:end);
        prevModR{iTrials} = prevModR{iTrials}(:,1:gaussShift:end);
    end
    
    %determine timepoint of reward given
    waterR{iTrials} = false(frames, round(sRate)*2);
    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = round((preStimDur + water(iTrials)) * sRate); %timepoint in frames when reward was given
        waterR{iTrials}(:, 1: size(timeR,2) - waterOn + 1) = timeR(:, waterOn:end);
    end
    
    if gaussShift > 1
        waterR{iTrials} = waterR{iTrials}(:,1:gaussShift:end); % subsample regressor
    end
    
    %% lever grab / release
    lGrabR{iTrials} = false(frames, length(motorIdx));
    lGrabRelR{iTrials} = false(frames, length(motorIdx));
    
    rGrabR{iTrials} = false(frames, length(motorIdx));
    rGrabRelR{iTrials} = false(frames, length(motorIdx));
    
    [grabL,grabRelL] = checkLevergrab(tapDur,postStimDur,levGrabL{iTrials},levReleaseL{iTrials},1/sRate);
    grabL(grabL < -preStimDur) = []; %ensure there are no grab events too early in the trial
    if ~isempty(grabL);grabRelL(grabRelL<grabL(1)) = []; end %make sure there are no release events before first grab
    
    [grabR,grabRelR] = checkLevergrab(tapDur,postStimDur,levGrabR{iTrials},levReleaseR{iTrials},1/sRate);
    grabR(grabR < -preStimDur) = [];
    if ~isempty(grabR); grabRelR(grabRelR<grabR(1)) = []; end
    
    for iRegs = 0 : length(motorIdx)-1
        
        shiftGrabL = grabL - (mPreTime - (iRegs * 1/sRate));
        shiftGrabR = grabR - (mPreTime - (iRegs * 1/sRate));
        
        shiftGrabRelL = grabRelL - (mPreTime - (iRegs * 1/sRate));
        shiftGrabRelR = grabRelR - (mPreTime - (iRegs * 1/sRate));
        
        if iRegs == 0 %first regressor, find first grab / tap
            
            lGrabR{iTrials}(find(histcounts(shiftGrabL,-preStimDur:1/sRate:postStimDur),1),iRegs+1) = true;
            rGrabR{iTrials}(find(histcounts(shiftGrabR,-preStimDur:1/sRate:postStimDur),1),iRegs+1) = true;
            
        else
            
            regOut = leverEvents(iRegs, grabRelL, shiftGrabL, find(lGrabR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, 'pre'); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            lGrabR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabRelR, shiftGrabR, find(rGrabR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, 'pre'); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            rGrabR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabL, shiftGrabRelL, find(lGrabRelR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, ''); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            lGrabRelR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabR, shiftGrabRelR, find(rGrabRelR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, ''); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            rGrabRelR{iTrials}(regOut,iRegs+1) = true;           
            
        end
    end
        
    if gaussShift > 1
        % subsample regressors
        lGrabR{iTrials} = lGrabR{iTrials}(:,1:gaussShift:end);
        lGrabRelR{iTrials} = lGrabRelR{iTrials}(:,1:gaussShift:end);
        rGrabR{iTrials} = rGrabR{iTrials}(:,1:gaussShift:end);
        rGrabRelR{iTrials} = rGrabRelR{iTrials}(:,1:gaussShift:end);
    end
    
    %% pupil / whisk / nose regressors
    bhvFrameRate = round(1/mean(diff(pTime{bTrials(iTrials)})));
    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab(iTrials) - preStimDur);
    trialTime = pTime{bTrials(iTrials)} - trialOn;    
    idx = trialTime < trialDur; %don't use late frames
    trialTime = trialTime(idx);
    timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial

    if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
        addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
        trialTime = [trialTime' addTime];
    end
    
    fastPupilR{iTrials} = Behavior_vidResamp(fPupil{bTrials(iTrials)}(idx), trialTime, sRate);
    fastPupilR{iTrials} = fastPupilR{iTrials}(end - frames + 1 : end);
        
    slowPupilR{iTrials} = Behavior_vidResamp(sPupil{bTrials(iTrials)}(idx), trialTime, sRate);
    slowPupilR{iTrials} = slowPupilR{iTrials}(end - frames + 1 : end);
        
    whiskR{iTrials} = Behavior_vidResamp(whisker{bTrials(iTrials)}(idx), trialTime, sRate);
    whiskR{iTrials} = smooth(whiskR{iTrials}(end - frames + 1 : end), 'rlowess');
        
    noseR{iTrials} = Behavior_vidResamp(nose{bTrials(iTrials)}(idx), trialTime, sRate);
    noseR{iTrials} = smooth(noseR{iTrials}(end - frames + 1 : end), 'rlowess');
            
    %% give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
    
    % get piezo, aligned to handle grab
    piezoR{iTrials} = data.analog(piezoLine,aCenter - floor((preStimDur + stimTime(1)) * data.analogFreq) + 1 : aCenter + floor((postStimDur - stimTime(1)) * data.analogFreq), iTrials);
end

%% rebuild whisker motion to get proper design matrix
wTrace = double(cat(1,whiskR{:}));
wTrace = (wTrace - prctile(wTrace,1))./ nanstd(wTrace); %minimum values are at 0, signal in standard deviation units
aWhisk = wTrace; %keep analog trace
wTrace = wTrace > 2; %take activity above 4 SDUs as indicator for whisking
wTrace = diff([0; wTrace]) == 1; %find event onsets
wTrace = reshape(wTrace,[],trialCnt);

for iTrials = 1:trialCnt
    
    trace = logical(histcounts(find(wTrace(:,iTrials)),0:bhvFrameRate/round(sRate):(bhvFrameRate/round(sRate))*frames))'; %resample to imaging frame rate. This is the zero lag regressor.    
    licks = round(([lickL{iTrials} lickR{iTrials}] + preStimDur) * round(sRate));  %lick times in frames. Dont check for whisking there because licking may move the snout around as well.
    if ~isempty(licks)
        lickZone = reshape(bsxfun(@plus, licks', -noWhisk : noWhisk),1,[]); %done use whisks that are too close to licks.
        trace(lickZone) = false;
    end
    
    % create full design matrix
    cIdx = bsxfun(@plus,find(trace),motorIdx);
    cIdx(cIdx < 1) = 0;
    cIdx(cIdx > frames) = frames;
    cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(motorIdx)-1));
    cIdx(cIdx < 1) = frames;
    cIdx(cIdx > (frames * length(motorIdx))) = frames * length(motorIdx);
    
    whiskR{iTrials} = false(frames, length(motorIdx));
    whiskR{iTrials}(cIdx(:)) = true;
    whiskR{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
    
    if gaussShift > 1
        whiskR{iTrials} = whiskR{iTrials}(:,1:gaussShift:end);
    end
end
whiskR = [single(aWhisk) cat(1,whiskR{:})]; %combine trials
clear wTrace idx cIdx aWhisk

%% rebuild nose motion to get proper design matrix
nTrace = double(cat(1,noseR{:}));
nTrace = (nTrace - prctile(nTrace,1))./ nanstd(nTrace); %minimum values are at 0, signal in standard deviation units
aSniff = nTrace; %keep analog trace
nTrace = nTrace > 2; %take activity above 4 SDUs as indicator for sniffing
nTrace = diff([0; nTrace]) == 1; %find event onsets
nTrace = reshape(nTrace,[],trialCnt);

for iTrials = 1:trialCnt
    
    trace = logical(histcounts(find(nTrace(:,iTrials)),0:bhvFrameRate/round(sRate):(bhvFrameRate/round(sRate))*frames))'; %resample to imaging frame rate. This is the zero lag regressor.    
    licks = round(([lickL{iTrials} lickR{iTrials}] + preStimDur) * round(sRate));  %lick times in frames. Dont check for sniffing there because licking may move the snout around as well.
    if ~isempty(licks)
        lickZone = reshape(bsxfun(@plus, licks', -noWhisk : noWhisk),1,[]); %done use whisks that are too close to licks.
        trace(lickZone) = false;
    end
    
    % create full design matrix
    cIdx = bsxfun(@plus,find(trace),motorIdx);
    cIdx(cIdx < 1) = 0;
    cIdx(cIdx > frames) = frames;
    cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(motorIdx)-1));
    cIdx(cIdx < 1) = frames;
    cIdx(cIdx > (frames * length(motorIdx))) = frames * length(motorIdx);
    
    noseR{iTrials} = false(frames, length(motorIdx));
    noseR{iTrials}(cIdx(:)) = true;
    noseR{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
    
    if gaussShift > 1
        noseR{iTrials} = noseR{iTrials}(:,1:gaussShift:end);
    end
end
noseR = [single(aSniff) cat(1,noseR{:})]; %combine trials
clear nTrace idx cIdx aSniff

%% rebuild piezo sensor to get proper design matrix
pTrace = zscore(double(cat(2,piezoR{:})));
pTrace = reshape(pTrace,[],trialCnt);
aPiezo = abs([zeros(1,trialCnt);diff(resample(pTrace,round(sRate),1000))]);
aPiezo = aPiezo(:) - prctile(aPiezo(:),1);

for iTrials = 1:trialCnt
    
    trace = smooth(pTrace(:,iTrials),round(sRate),'lowess');
    trace = trace > 0.5; %threshold normalized sensor data
    trace = diff([0;trace]) == 1; %find event onsets
    trace = logical(histcounts(find(trace), 0:1000/round(sRate):(1000/round(sRate))*frames))'; %resample to imaging frame rate. This is the zero lag regressor.
    trace(1:5) = false; %first frames should not be onset events to avoid misinterpration of onging sensor read on first datapoint.
    cIdx = bsxfun(@plus,find(trace),motorIdx);
    cIdx(cIdx < 1) = 0;
    cIdx(cIdx > frames) = frames;
    cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(motorIdx)-1));
    cIdx(cIdx < 1) = frames;
    cIdx(cIdx > (frames * length(motorIdx))) = frames * length(motorIdx);
    
    piezoR{iTrials} = false(frames, length(motorIdx));
    piezoR{iTrials}(cIdx(:)) = true;
    piezoR{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
    piezoR{iTrials}(end,2:end) = piezoR{iTrials}(end-1,1:end-1); %replace with shifted version of previous timepoint
    
    if gaussShift > 1
        piezoR{iTrials} = piezoR{iTrials}(:,1:gaussShift:end);
    end
end
piezoR = [single(aPiezo) cat(1,piezoR{:})]; %combine trials
clear pTrace cIdx aPiezo

%% re-align behavioral video data and neural data to lever grab instead of stimulus onset
V1 = reshape(V1,205,[],bhvDimCnt); %get to trial format
vidR = V1(:,bTrials,:); clear V1 %get correct trials from behavioral video data.

% re-align video data
temp1 = NaN(size(data.dFOF,1),frames,trialCnt);
temp2 = NaN(frames,trialCnt,bhvDimCnt);
for x = 1 : size(vidR,2)
    try
        temp1(:,:,x) = data.dFOF(:,(94 - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (94 - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        temp2(:,x,:) = vidR((94 - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (94 - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
    catch
        fprintf(1,'Could not align trial %d. Relative stim time: %fs\n', x, stimTime(x));
    end
end
Vc = reshape(temp1,size(data.dFOF,1),[]); clear temp1
vidR = reshape(temp2,[],bhvDimCnt); clear temp2

%% reshape regressors, make design matrix and indices for regressors that are used for the model
timeR = repmat(logical(diag(ones(1,frames))),trialCnt,1); %time regressor

lGrabR = cat(1,lGrabR{:});
lGrabRelR = cat(1,lGrabRelR{:});
rGrabR = cat(1,rGrabR{:});
rGrabRelR = cat(1,rGrabRelR{:});

lLickR = cat(1,lLickR{:});
rLickR = cat(1,rLickR{:});
leverInR = cat(1,leverInR{:});

lVisStimR = cat(1,lVisStimR{:});
rVisStimR = cat(1,rVisStimR{:});
lAudStimR = cat(1,lAudStimR{:});
rAudStimR = cat(1,rAudStimR{:});

visRewardR = cat(1,visRewardR{:});
audRewardR = cat(1,audRewardR{:});
prevRewardR = cat(1,prevRewardR{:});

visChoiceR = cat(1,visChoiceR{:});
audChoiceR = cat(1,audChoiceR{:});
prevChoiceR = cat(1,prevChoiceR{:});

prevModR = cat(1,prevModR{:});
waterR = cat(1,waterR{:});

fastPupilR = cat(1,fastPupilR{:});
fastPupilR(~isnan(fastPupilR(:,1)),:) = zscore(fastPupilR(~isnan(fastPupilR(:,1)),:));

slowPupilR = cat(1,slowPupilR{:});
slowPupilR(~isnan(slowPupilR(:,1)),:) = zscore(slowPupilR(~isnan(slowPupilR(:,1)),:));

if ~any(isnan(fastPupilR(:,1)) == isnan(vidR(:,1)))
    error('Pupil and bhv video do not agree with one another. Maybe something wrong with timing ?')
end

%% get video motion and create trial-aligned motion regressors
% orthogonalize video against all other regressors before computing motion
orthVidR = NaN(size(vidR));
smallR = [timeR visChoiceR audChoiceR visRewardR audRewardR lGrabR lGrabRelR rGrabR rGrabRelR lLickR rLickR leverInR lVisStimR rVisStimR lAudStimR rAudStimR ...
    prevRewardR prevChoiceR prevModR waterR piezoR whiskR noseR fastPupilR slowPupilR vidR];
trialIdx = ~isnan(mean(smallR,2));
smallR = smallR(trialIdx,:);
smallR(:,nansum(abs(smallR)) < 10) = []; %reject regressors that are too sparse

[Q, redQRR] = qr(bsxfun(@rdivide,smallR,sqrt(sum(smallR.^2))),0); %orthogonalize video against small design matrix
orthVidR(trialIdx,:) = Q(:,end-size(vidR,2)+1:end); clear smallR

vidMove = abs(diff(orthVidR));%change in video data. This is our proxy for animal motion.
vidMove = [vidMove(1,:);vidMove];
vidMove(frames+1 : frames : size(vidMove,1)-frames+1, :) = vidMove(frames+2 : frames : size(vidMove,1)-frames+2, :); %remove trial transition transient
moveR = vidMove; %keep analog traces
vidMove = nanmean(reshape(vidMove,frames,[],size(vidR,2)),3); %reshape to trials and average over dimensions

for x = 1: length(trialSegments)-1 %this should be 4 segments in total (Baseline, Handle, Stim and Waiting period)
   trialMove(x,:) = mean(vidMove(trialSegments{x},:));  %compute average motion change in each segment for each trial
end

%compute average motion change in each segment for each trial
for iTrials = 1:size(trialMove,2)
    BaselineMoveR{iTrials} = diag(repmat(trialMove(1,iTrials),frames,1));
    HandleMoveR{iTrials} = diag(repmat(trialMove(2,iTrials),frames,1));
    StimulusMoveR{iTrials} = diag(repmat(trialMove(3,iTrials),frames,1));
    WaitMoveR{iTrials} = diag(repmat(trialMove(4,iTrials),frames,1));
end

BaselineMoveR = cat(1,BaselineMoveR{:});
HandleMoveR = cat(1,HandleMoveR{:});
StimulusMoveR = cat(1,StimulusMoveR{:});
WaitMoveR = cat(1,WaitMoveR{:});

%% create full design matrix
fullR = [timeR visChoiceR audChoiceR visRewardR audRewardR lGrabR lGrabRelR rGrabR rGrabRelR lLickR rLickR leverInR lVisStimR rVisStimR lAudStimR rAudStimR ...
    prevRewardR prevChoiceR prevModR waterR piezoR whiskR noseR fastPupilR slowPupilR BaselineMoveR HandleMoveR StimulusMoveR WaitMoveR moveR vidR];

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
recLabels = {
    'time' 'visChoice' 'audChoice' 'visReward' 'audReward' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'leverIn' 'lVisStim' 'rVisStim' ...
    'lAudStim' 'rAudStim' 'prevReward' 'prevChoice' 'prevMod' 'water' 'piezo' 'whisk' 'nose' 'fastPupil' 'slowPupil' 'BaselineMove' 'HandleMove' 'StimulusMove' 'WaitMove' 'Move' 'bhvVideo'};

%index to reconstruct different response kernels
recIdx = [
    ones(1,size(timeR,2))*find(ismember(recLabels,'time')) ...
    ones(1,size(visChoiceR,2))*find(ismember(recLabels,'visChoice')) ...
    ones(1,size(audChoiceR,2))*find(ismember(recLabels,'audChoice')) ...
    ones(1,size(visRewardR,2))*find(ismember(recLabels,'visReward')) ...
    ones(1,size(audRewardR,2))*find(ismember(recLabels,'audReward')) ...
    ones(1,size(lGrabR,2))*find(ismember(recLabels,'lGrab')) ...
    ones(1,size(lGrabRelR,2))*find(ismember(recLabels,'lGrabRel')) ...
    ones(1,size(rGrabR,2))*find(ismember(recLabels,'rGrab')) ...
    ones(1,size(rGrabRelR,2))*find(ismember(recLabels,'rGrabRel')) ...
    ones(1,size(lLickR,2))*find(ismember(recLabels,'lLick')) ...
    ones(1,size(rLickR,2))*find(ismember(recLabels,'rLick')) ...
    ones(1,size(leverInR,2))*find(ismember(recLabels,'leverIn')) ...
    ones(1,size(lVisStimR,2))*find(ismember(recLabels,'lVisStim')) ...
    ones(1,size(rVisStimR,2))*find(ismember(recLabels,'rVisStim')) ...
    ones(1,size(lAudStimR,2))*find(ismember(recLabels,'lAudStim')) ...
    ones(1,size(rAudStimR,2))*find(ismember(recLabels,'rAudStim')) ...
    ones(1,size(prevRewardR,2))*find(ismember(recLabels,'prevReward')) ...
    ones(1,size(prevChoiceR,2))*find(ismember(recLabels,'prevChoice')) ...
    ones(1,size(prevModR,2))*find(ismember(recLabels,'prevMod')) ...
    ones(1,size(waterR,2))*find(ismember(recLabels,'water')) ...
    ones(1,size(piezoR,2))*find(ismember(recLabels,'piezo')) ...
    ones(1,size(whiskR,2))*find(ismember(recLabels,'whisk')) ...
    ones(1,size(noseR,2))*find(ismember(recLabels,'nose')) ...
    ones(1,size(fastPupilR,2))*find(ismember(recLabels,'fastPupil')) ...
    ones(1,size(slowPupilR,2))*find(ismember(recLabels,'slowPupil')) ...
    ones(1,size(BaselineMoveR,2))*find(ismember(recLabels,'BaselineMove')) ...
    ones(1,size(HandleMoveR,2))*find(ismember(recLabels,'HandleMove')) ...
    ones(1,size(StimulusMoveR,2))*find(ismember(recLabels,'StimulusMove')) ...
    ones(1,size(WaitMoveR,2))*find(ismember(recLabels,'WaitMove')) ...
    ones(1,size(moveR,2))*find(ismember(recLabels,'Move')) ...
    ones(1,size(vidR,2))*find(ismember(recLabels,'bhvVideo'))];

trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fprintf(1, 'Rejected %d/%d trials for NaN entries in regressors\n', sum(trialIdx)/frames,trialCnt);
fullR(trialIdx,:) = []; %clear bad trials

idx = nansum(abs(fullR)) < 10; %reject regressors that are too sparse
fullR(:,idx) = []; %clear empty regressors
fprintf(1, 'Rejected %d/%d empty regressors\n', sum(idx),length(idx));

%% run QR and check for rank-defficiency
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize design matrix
figure; plot(abs(diag(fullQRR))); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    error('Design matrix is rank-defficient')
end

%% save modified Vc
Vc(:,trialIdx) = []; %clear bad trials
save([cPath 'interpVc.mat'], 'Vc', 'frames');

%% apply gaussian filter to design matrix if using sub-sampling
if gaussShift > 1
    [a,b] = size(fullR);
    
    % find non-continous regressors (contain values different from -1, 0 or 1)
    temp = false(size(fullR));
    temp(fullR(:) ~= 0 & fullR(:) ~= 1 & fullR(:) ~= -1 & ~isnan(fullR(:))) = true;
    regIdx = nanmean(temp) == 0; %index for non-continous regressors
    
    % do gaussian convolution. perform trialwise to avoid overlap across trials.
    trialCnt = a/frames;
    fullR = reshape(fullR,frames,trialCnt,b);
    for iTrials = 1:trialCnt
        fullR(:,iTrials,regIdx) = smoothCol(squeeze(fullR(:,iTrials,regIdx)),gaussShift*2,'gauss');
    end
    fullR = reshape(fullR,a,b);
end

%% clear individual regressors
clear stimR lGrabR lGrabRelR rGrabR rGrabRelR waterR lLickR rLickR leverInR ...
    lVisStimR rVisStimR lAudStimR rAudStimR visRewardR audRewardR prevRewardR visChoiceR audChoiceR ...
    prevChoiceR prevModR fastPupilR moveR BaselineMoveR HandleMoveR StimulusMoveR WaitMoveR piezoR whiskR noseR

%% run ridge regression in low-D, with and without ortogonalized video
%run model with non-orthogonalized video. Zero-mean without intercept.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'orgdimBeta.mat'], 'dimBeta', 'ridgeVals')
save([cPath filesep 'orgregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','fullQRR','-v7.3');
Behavior_betaRebuild(cPath, 'org'); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

%run model with orthogonalized video. Zero-mean without intercept.
vidIdx = ismember(recIdx(~idx), find(ismember(recLabels,{'bhvVideo'}))); % index for video regressors
fullR(:,vidIdx) = orthVidR(~trialIdx,:); %orthogonalize video regressors
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals')
save([cPath filesep 'regData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','redQRR','-v7.3');
Behavior_betaRebuild(cPath); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data


%% Re-run model with baseline correction and intercept
% compute baseline-corrected Vc
Vc = reshape(Vc, size(data.dFOF,1), frames, []);
temp = squeeze(mean(mean(Vc(:,1:sRate/2,:),3),2));
Vc = reshape(Vc, size(data.dFOF,1), []);
Vc = bsxfun(@minus, Vc, temp);

%run model with orthogonalized video. Baseline-corrected, with intercept.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, false); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty non zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'offsetdimBeta.mat'], 'dimBeta', 'ridgeVals')

fullR = [ones(size(fullR,1), 1) fullR]; %add intercept to design matrix
recIdx = [1 recIdx+1]; %add zero for intercept
idx = [false idx]; %add zero for intercept
recLabels = ['offset' recLabels];
save([cPath filesep 'offsetregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','redQRR','-v7.3');

%% run video only model. Zero-mean without intercept.
fullR = vidR(~trialIdx,:);
idx = false(1,size(vidR,2));
recLabels = {'bhvVideo'};
recIdx = ones(1,size(vidR,2));
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for video-only zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'vidOnlydimBeta.mat'], 'dimBeta', 'ridgeVals')
save([cPath filesep 'vidOnlyregData.mat'], 'fullR', 'idx', 'trialIdx', 'recIdx', 'recLabels','gaussShift','-v7.3');
Behavior_betaRebuild(cPath, 'vidOnly'); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data


%% nested functions
function regOut = leverEvents(cReg, times, shiftTimes, lastTimes, timeFrame, preTimeShift, dMode) %code to compute events for current lever regressor

regOut = [];
if cReg == 0 %first regressor, find first grab
    regOut = find(histcounts(shiftTimes,timeFrame),1);
    
elseif cReg < preTimeShift && strcmpi(dMode,'pre') %negative shift regressors - only for grab or tap
    regOut = lastTimes + 1; %find last event and shift one forward
    if isempty(regOut)
        regOut = find(histcounts(shiftTimes,timeFrame),1); %check if new event can be found if none is present so far
    end
    
elseif cReg == preTimeShift %this is the zero-lag regressor. use all available events
    regOut = find(histcounts(shiftTimes,timeFrame)); %times and shifttimes should be the same here.
    
elseif cReg > preTimeShift %for positive shift, use the last regressor but eliminate event if there is overlap with zero-lag regressor.
    regOut = lastTimes + 1; %find last event and shift one forward
    if ~isempty(regOut)
        regOut(ismember(regOut, find(histcounts(times,timeFrame)))) = []; %remove events that overlap with the given zero-lag regressor. Can also be from a different event type like handle release.
    end
end



function [grabOn,grabRel,tapOn,tapRel] = checkLevergrab(tapDur,postStimDur,grabs,release,minTime)

grabOn = [];
grabRel = [];
tapOn = [];
tapRel = [];

if ~isempty(grabs)
    Cnt = 0;
    grabCnt = 0;
    tapCnt = 0;
    
    while true
        Cnt = Cnt +1;
        if length(grabs) < Cnt %no more grabs
            break;
        end
        
        idx = find(release > grabs(Cnt),1);
        if ~isempty(idx) %if there is a release that follows current grab
            cGrab = release(idx) - grabs(Cnt); %duration of current grab
        else
            cGrab = postStimDur - grabs(Cnt);
        end
        
        %% more grabs available - check start time of next grab and merge if they are very close in time
        if length(grabs) > Cnt && ~isempty(idx)
            while (grabs(Cnt+1) - release(idx)) <= minTime %time diff between grabs is less than minimum duration. Merge with next grab.
                
                release(idx) = []; %delete current release
                grabs(Cnt+1) = []; %delete next grab
                
                idx = find(release > grabs(Cnt),1); %next release
                if ~isempty(idx) %if there is a release that follows current grab
                    cGrab = release(idx) - grabs(Cnt); %duration of current grab
                else
                    cGrab = postStimDur - grabs(Cnt);
                end
                
                if length(grabs) <= Cnt %no more grabs
                    break
                end
                
            end
        end
        
        %% check if current grab is grab or tap
        if cGrab <= tapDur
            tapCnt = tapCnt + 1;
            tapOn(tapCnt) = grabs(Cnt);
            if ~isempty(idx)
                tapRel(tapCnt) = idx;
            end
        else
            grabCnt = grabCnt + 1;
            grabOn(grabCnt) = grabs(Cnt);
            if ~isempty(idx)
                grabRel(grabCnt) = release(idx);
            end
        end
        
        if isempty(idx) || length(grabs) <= Cnt %no more grabs/releases
            break;
        end
    end
end