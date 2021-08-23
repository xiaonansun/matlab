function rateDisc_choiceModelDLC(cPath,Animal,Rec,dType)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

if ~exist('dType', 'var') || isempty(dType)
    dType = 'Widefield'; %default is widefield data
end

if strcmpi(dType,'twoP')
    piezoLine = 5;     % channel in the analog data that contains data from piezo sensor
    stimLine = 4;      % channel in the analog data that contains stimulus trigger.
    
elseif strcmpi(dType,'Widefield')
    piezoLine = 2;     % channel in the analog data that contains data from piezo sensor
    stimLine = 6;      % channel in the analog data that contains stimulus trigger.
end
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
if ispc
    sPath = ['\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
else
    sPath = ['/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
end
disp(cPath);

%% general variables
load([cPath 'opts2.mat'], 'opts'); % get some options from imaging data
if opts.frameRate == 15
    vcFile = 'Vc.mat';
elseif opts.frameRate == 30 %if data was acquired at 30Hz, there should be a resampled Vc (to 15Hz)
    vcFile = 'rsVc.mat';
    opts.frameRate = 15;
end

sRate = opts.frameRate; % Sampling rate of imaging in Hz
preStimDur = floor(2 * sRate) / sRate; % Duration of trial before lever grab in seconds
postStimDur = floor(3 *sRate) / sRate; % Duration of trial after lever grab onset in seconds
frames = round((preStimDur + postStimDur) * sRate); %nr of frames per trial
trialDur = (frames * (1/sRate)); %duration of trial in seconds

%other variables
opts.mPreTime = ceil(0.5 * sRate);  % precede motor events to capture preparatory activity in frames
opts.mPostTime = ceil(2 * sRate);   % follow motor events for mPostStim in frames
fsPreTime = ceil(1 * sRate);   % preceed first stim event for fsPreTime in frames
fsPostTime = ceil(2 * sRate);   % follow first stim event for fsPostStim in frames
fsPostTime = fsPreTime + fsPostTime; %total time for stimulus kernel
sPostTime = ceil(1 * sRate);   % follow stim events for sPostStim in frames
motorIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post motor action
tapDur = 0.1;      % minimum time of lever contact, required to count as a proper grab.
leverMoveDur = 0.3; %duration of lever movement. this is used to orthogonalize video against lever movement.
leverMoveDur = ceil(leverMoveDur * sRate); %convert to frames
ridgeFolds = 10;    %number of cross-validations for motor/task models
shVal = sRate * opts.preStim  + 1; %expected position of stimulus onset in the imaging data (s).
maxStimShift = 1 * sRate; % maximal stimulus onset after handle grab. (default is 1s - this means that stimulus onset should not be more than 1s after handle grab. Cognitive regressors will have up to 1s of baseline because of stim jitter.)
bhvDimCnt = 200;    % number of dimensions from behavioral videos that are used as regressors.
gaussShift = 1;     % inter-frame interval between regressors. Will use only every 'gaussShift' regressor and convolve with gaussian of according FHWM to reduce total number of used regressors.
[~, motorLabels] = rateDiscRecordings; %get motor labels for motor-only model
opMotorLabels = {'lLick' 'rLick' 'lGrab' 'rGrab'}; %operant motor regressors
dims = 200; %number of dims in Vc that are used in the model
opts.framesPerTrial = frames;
dlcThresh = 1; %threshold for event detection for deeplabcut labels in SDUs
lickThresh = 0.5; %threshold for lick detection for deeplabcut tongue labels (SDUs)
pupilThresh = 1; %threshold for pupil dilation/constriction in SDUs
pupilShift = 0.75; %shift back pupil measure to compensate for muscle delay

%% load data
bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data
SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds

if strcmpi(dType,'Widefield')
    if exist([cPath vcFile],'file') ~= 2 %check if data file exists and get from server otherwise
        copyfile([sPath vcFile],[cPath vcFile]);
        bhvFile = dir([sPath filesep Animal '_' Paradigm '*.mat']);
        copyfile([sPath bhvFile.name],[cPath bhvFile.name]);
    end
    load([cPath vcFile], 'Vc', 'U', 'trials', 'bTrials')
    Vc = Vc(1:dims,:,:);
    U = U(:,:,1:dims);
    
    % ensure there are not too many trials in Vc
    ind = trials > SessionData.nTrials;
    trials(ind) = [];
    bTrials(ind) = [];
    Vc(:,:,ind) = [];

elseif strcmpi(dType,'twoP')
    
    load([cPath 'data'], 'data'); %load 2p data
    % ensure there are not too many trials in the dataset
    bTrials = data.trialNumbers;
    trials = bTrials;
    bTrials(~ismember(data.trialNumbers,data.bhvTrials)) = []; %don't use trials that have problems with trial onset times
    bTrials(SessionData.DidNotChoose(bTrials) | SessionData.DidNotLever(bTrials) | ~SessionData.Assisted(bTrials)) = []; %don't use unperformed/assisted trials
    
    data.dFOF(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    data.DS(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    data.analog(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    
    Vc = data.dFOF; %Vc is now neurons x frames x trials
    dims = size(data.dFOF,1); %dims is now # of neurons instead
end
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset

%equalize L/R choices
useIdx = ~isnan(bhv.ResponseSide); %only use performed trials
choiceIdx = rateDisc_equalizeTrials(useIdx, bhv.ResponseSide == 1, bhv.Rewarded, inf, true); %equalize correct L/R choices
trials = trials(choiceIdx);
bTrials = bTrials(choiceIdx);
Vc = Vc(:,:,choiceIdx);

bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset   
trialCnt = length(bTrials);

%% load behavior data
load([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'vidV'); %load behavior video data
V1 = vidV(:,1:bhvDimCnt); %behavioral video regressors
load([cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'vidV'); %load abs motion video data
V2 = vidV(:,1:bhvDimCnt); % motion regressors

% load dlc data
load([cPath 'dlcData.mat'], 'dlcData', 'dlcVars', 'dlcRate'); %load deeplabcut labels
for iCams = 1 : 2
    dlcData{iCams} = reshape(cat(1,dlcData{iCams}{:,:}),[],length(dlcVars{iCams})); %convert to matrix
    dlcVars{iCams} = cellfun(@(x) ['cam' num2str(iCams) '_' x], dlcVars{iCams},'UniformOutput',false); %add camera to variable names
end
%combine cameras
dlcData = cat(2,dlcData{:});
dlcVars = cat(2,dlcVars{:});

% check options that were used for dimensionality reduction and ensure that imaging and video data trials are equal length
load([cPath 'BehaviorVideo' filesep 'bhvOpts.mat'], 'bhvOpts'); %load abs motion video data
bhvRate = bhvOpts.targRate; %framerate of face camera
if bhvRate ~= dlcRate
    error('Sampling rate of video and DLC dont match.')
end
if (bhvOpts.preStimDur + bhvOpts.postStimDur) > (opts.preStim + opts.postStim) %if behavioral video trials are longer than imaging data
    V1 = reshape(V1, [], SessionData.nTrials, bhvDimCnt);
    V2 = reshape(V2, [], SessionData.nTrials, bhvDimCnt);
    dlcData = reshape(dlcData, [], SessionData.nTrials, length(dlcVars));
    
    if bhvOpts.preStimDur > opts.preStim
        frameDiff = ceil((bhvOpts.preStimDur - opts.preStim) * bhvRate); %overhead in behavioral frames that can be removed.
        V1 = V1(frameDiff+1:end, :, :); %cut data to size
        V2 = V2(frameDiff+1:end, :, :);
        dlcData = dlcData(frameDiff+1:end, :, :);
    end
    if bhvOpts.postStimDur > opts.postStim
        frameDiff = ceil((bhvOpts.postStimDur - opts.postStim) * bhvRate); %overhead in behavioral frames that can be removed.
        V1 = V1(1 : end - frameDiff, :, :); %cut data to size
        V2 = V2(1 : end - frameDiff, :, :);
        dlcData = dlcData(1 : end - frameDiff, :, :);
    end
    V1 = reshape(V1, [], bhvDimCnt);
    V2 = reshape(V2, [], bhvDimCnt);
    dlcData = reshape(dlcData, [], length(dlcVars));
end

load([cPath 'BehaviorVideo' filesep 'FilteredPupil.mat'], 'pTime', 'fPupil', 'whisker', 'faceM', 'bodyM', 'nose', 'bTime'); %load pupil data
%check if timestamps from pupil data are shifted against bhv data
timeCheck1 = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
timeCheck2 = (SessionData.TrialStartTime(1)) - (bTime{1}(1)); %time difference between first acquired frame and onset of first trial
if (timeCheck1 > 3590 && timeCheck1 < 3610) && (timeCheck2 > 3590 && timeCheck2 < 3610) %timeshift by one hour (+- 10seconds)
    warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
    for iTrials = 1 : length(pTime)
        pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
        bTime{iTrials} = bTime{iTrials} + 3600; %add one hour
    end
elseif timeCheck1 > 30 || timeCheck1 < -30 || timeCheck2 > 30 || timeCheck2 < -30
    error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
end

if any(bTrials > length(pTime))
    warning(['There are insufficient trials in the pupil data. Rejected the last ' num2str(sum(bTrials > length(pTime))) ' trial(s)']);
    bTrials(bTrials > length(pTime)) = [];
    trialCnt = length(bTrials);
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
handleSounds = cell(1,trialCnt);

tacStimL = cell(1,trialCnt);
tacStimR = cell(1,trialCnt);
audStimL = cell(1,trialCnt);
audStimR = cell(1,trialCnt);

for iTrials = 1:trialCnt
    
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    try
        stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
        handleSounds{iTrials} = leverTimes(1:2:end) - stimGrab(iTrials); %track indicator sound when animal is grabing both handles
        stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab(iTrials); %time of stimulus onset - measured from soundcard
        stimEndTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.DecisionWait(1) - stimGrab(iTrials); %end of stimulus period, relative to handle grab
    catch
        stimTime(iTrials) = NaN;
        stimEndTime(iTrials) = NaN;
        stimGrab(iTrials) = 0;
    end
    
    %check for spout motion
    if isfield(bhv.RawEvents.Trial{iTrials}.States,'MoveSpout')
        spoutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - stimGrab(iTrials);
        
        %also get time when the other spout was moved out at
        if bhv.Rewarded(iTrials)
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
        else
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.HardPunish(1) - stimGrab(iTrials);
        end
    else
        spoutTime(iTrials) = NaN;
        spoutOutTime(iTrials) = NaN;
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
    
    % get stimulus events times
    audStimL{iTrials} = bhv.stimEvents{iTrials}{1} + stimTime(iTrials);
    audStimR{iTrials} = bhv.stimEvents{iTrials}{2} + stimTime(iTrials);
    tacStimL{iTrials} = bhv.stimEvents{iTrials}{5} + stimTime(iTrials);
    tacStimR{iTrials} = bhv.stimEvents{iTrials}{6} + stimTime(iTrials);
    
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
maxSpoutRegs = length(min(round((preStimDur + spoutTime) * sRate)) : frames); %maximal number of required spout regressors

%% build regressors - create design matrix based on event times
%basic time regressors
timeR = logical(diag(ones(1,frames)));

lGrabR = cell(1,trialCnt);
rGrabR = cell(1,trialCnt);
lLickR = cell(1,trialCnt);
rLickR = cell(1,trialCnt);
leverInR = cell(1,trialCnt);

lfirstTacStimR = cell(1,trialCnt);
rfirstTacStimR = cell(1,trialCnt);
lfirstAudStimR = cell(1,trialCnt);
rfirstAudStimR = cell(1,trialCnt);

lTacStimR = cell(1,trialCnt);
rTacStimR = cell(1,trialCnt);
lAudStimR = cell(1,trialCnt);
rAudStimR = cell(1,trialCnt);

spoutR = cell(1,trialCnt);
spoutOutR = cell(1,trialCnt);

rewardR = cell(1,trialCnt);
lhandleChoiceR = cell(1,trialCnt);
rhandleChoiceR = cell(1,trialCnt);
lstimChoiceR = cell(1,trialCnt);
rstimChoiceR = cell(1,trialCnt);
lresponseChoiceR = cell(1,trialCnt);
rresponseChoiceR = cell(1,trialCnt);

prevRewardR = cell(1,trialCnt);
prevChoiceR = cell(1,trialCnt);

waterR = cell(1,trialCnt);
pupilR = cell(1,trialCnt);

whiskR = cell(1,trialCnt);
noseR = cell(1,trialCnt);
piezoR = cell(1,trialCnt);
bodyR = cell(1,trialCnt);

handleSoundR = cell(1,trialCnt);

%%
tic
for iTrials = 1:trialCnt
    %% first tactile/auditory stimuli
    lfirstAudStimR{iTrials} = false(frames, fsPostTime);
    rfirstAudStimR{iTrials} = false(frames, fsPostTime);
    lfirstTacStimR{iTrials} = false(frames, fsPostTime);
    rfirstTacStimR{iTrials} = false(frames, fsPostTime);
    
    firstStim = NaN;
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 6 %auditory or mixed stimulus
        if ~isempty(audStimL{iTrials}(~isnan(audStimL{iTrials})))
            firstStim = round((audStimL{iTrials}(1) + preStimDur) * sRate) - fsPreTime;
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            lfirstAudStimR{iTrials}(:,1 : stimEnd - firstStim + 1)  = timeR(:, firstStim : stimEnd);
        end
        if ~isempty(audStimR{iTrials}(~isnan(audStimR{iTrials})))
            firstStim = round((audStimR{iTrials}(1) + preStimDur) * sRate) - fsPreTime;
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            rfirstAudStimR{iTrials}(:,1 : stimEnd - firstStim + 1) = timeR(:, firstStim : stimEnd);
        end
    end
    
    if bhv.StimType(iTrials) == 4 || bhv.StimType(iTrials) == 6 %tactile or mixed stimulus
        if ~isempty(tacStimL{iTrials}(~isnan(tacStimL{iTrials})))
            firstStim = round((tacStimL{iTrials}(1) + preStimDur) * sRate) - fsPreTime;
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            lfirstTacStimR{iTrials}(:,1 : stimEnd - firstStim + 1) = timeR(:, firstStim : stimEnd);
        end
        if ~isempty(tacStimR{iTrials}(~isnan(tacStimR{iTrials})))
            firstStim = round((tacStimR{iTrials}(1) + preStimDur) * sRate) - fsPreTime;
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            rfirstTacStimR{iTrials}(:,1 : stimEnd - firstStim + 1) = timeR(:, firstStim : stimEnd);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        lfirstTacStimR{iTrials} = lfirstTacStimR{iTrials}(:,1:gaussShift:end);
        rfirstTacStimR{iTrials} = rfirstTacStimR{iTrials}(:,1:gaussShift:end);
        lfirstAudStimR{iTrials} = lfirstAudStimR{iTrials}(:,1:gaussShift:end);
        rfirstAudStimR{iTrials} = rfirstAudStimR{iTrials}(:,1:gaussShift:end);
    end
    
    %% other tactile/auditory stimuli
    lAudStimR{iTrials} = false(frames, sPostTime);
    rAudStimR{iTrials} = false(frames, sPostTime);
    lTacStimR{iTrials} = false(frames, sPostTime);
    rTacStimR{iTrials} = false(frames, sPostTime);
    
    for iRegs = 0 : sPostTime - 1
        allStims = audStimL{iTrials}(2:end) + (iRegs * 1/sRate);
        lAudStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        allStims = audStimR{iTrials}(2:end) + (iRegs * 1/sRate);
        rAudStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        allStims = tacStimL{iTrials}(2:end) + (iRegs * 1/sRate);
        lTacStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        allStims = tacStimR{iTrials}(2:end) + (iRegs * 1/sRate);
        rTacStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    if gaussShift > 1
        % subsample regressors
        lTacStimR{iTrials} = lTacStimR{iTrials}(:,1:gaussShift:end);
        rTacStimR{iTrials} = rTacStimR{iTrials}(:,1:gaussShift:end);
        lAudStimR{iTrials} = lAudStimR{iTrials}(:,1:gaussShift:end);
        rAudStimR{iTrials} = rAudStimR{iTrials}(:,1:gaussShift:end);
    end
    
    %% spout regressors
    spoutIdx = round((preStimDur + spoutTime(iTrials)) * sRate) : round((preStimDur + postStimDur) * sRate); %index for which part of the trial should be covered by spout regressors
    spoutR{iTrials} = false(frames, maxSpoutRegs);
    if ~isnan(spoutTime(iTrials))
        spoutR{iTrials}(:, 1:length(spoutIdx)) = timeR(:, spoutIdx);
    end
    
    spoutOutR{iTrials} = false(frames, 3);
    spoutOut = round((preStimDur + spoutOutTime(iTrials)) * sRate); %time when opposing spout moved out again
    if ~isnan(spoutOut) && spoutOut < (frames + 1)
        cInd = spoutOut : spoutOut + 2; cInd(cInd > frames) = [];
        temp = diag(ones(1,3));
        spoutOutR{iTrials}(cInd, :) = temp(1:length(cInd),:);
    end
    
    if gaussShift > 1
        % subsample regressors
        spoutR{iTrials} = spoutR{iTrials}(:,1:gaussShift:end);
        spoutOutR{iTrials} = spoutOutR{iTrials}(:,1:gaussShift:end);
    end
    
    %% lick regressors
    lLickR{iTrials} = false(frames, length(motorIdx));
    rLickR{iTrials} = false(frames, length(motorIdx));
    
    for iRegs = 0 : length(motorIdx)-1
        licks = lickL{iTrials} - ((opts.mPreTime/sRate) - (iRegs * 1/sRate));
        lLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        licks = lickR{iTrials} - ((opts.mPreTime/sRate) - (iRegs * 1/sRate));
        rLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    if gaussShift > 1
        % subsample regressors
        lLickR{iTrials} = lLickR{iTrials}(:,1:gaussShift:end);
        rLickR{iTrials} = rLickR{iTrials}(:,1:gaussShift:end);
    end
    
    %% lever in
    leverInR{iTrials} = false(frames, leverMoveDur);
    leverShift = round((preStimDur + leverIn(iTrials))* sRate); %timepoint in frames when lever moved in, relative to lever grab
    
    if ~isnan(leverShift)
        if leverShift > 0 %lever moved in during the recorded trial
            leverInR{iTrials}(leverShift : leverShift + leverMoveDur -1, :) = diag(ones(1, leverMoveDur));
        elseif (leverShift + leverMoveDur) > 0  %lever was moving before data was recorded but still moving at trial onset
            leverInR{iTrials}(1 : leverMoveDur + leverShift, :) = [zeros(leverMoveDur + leverShift, abs(leverShift)) diag(ones(1, leverMoveDur + leverShift))];
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        leverInR{iTrials} = leverInR{iTrials}(:,1:gaussShift:end);
    end
    
    %% dual-handle indicator sound
    handleSoundR{iTrials} = false(frames, sPostTime);
    for iRegs = 0 : sPostTime - 1
        allStims = handleSounds{iTrials}(1:end) + (iRegs * 1/sRate);
        allStims = allStims(~isnan(allStims)) + 0.001;
        handleSoundR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    %% choice and reward - use four different episodes for choice: handle, stim, delay and response
    stimTemp = false(frames, frames + maxStimShift);
    stimShift = round(stimTime(iTrials) * sRate); %amount of stimshift. move diagonal on x-axis accordingly.
    if (stimShift > maxStimShift) || isnan(stimTime(iTrials))
        stimTemp = NaN(frames, frames + maxStimShift); %don't use trial if stim onset is too late
    else
        stimTemp(:, end - stimShift - frames + 1 : end - stimShift) = timeR;
    end
    
    rewardR{iTrials} = false(size(stimTemp));
    if bhv.Rewarded(iTrials) %rewarded
        rewardR{iTrials} = stimTemp; %trial was rewarded
    end
    
    % get L/R choices as binary design matrix. This includes 4 choice regressors before handle grab and subsequent regressors until stimulus onset
    handleReg = preStimDur*sRate; %regressor that is aligned with handle grab
    lhandleChoiceR{iTrials} = false(size(stimTemp,1),sRate);
    rhandleChoiceR{iTrials} = false(size(stimTemp,1),sRate);
    if stimShift > size(lhandleChoiceR{iTrials},2); stimShift = size(lhandleChoiceR{iTrials},2); end %make sure stimshift is not beyond size of design matrix
    try
        if bhv.ResponseSide(iTrials) == 1
            lhandleChoiceR{iTrials}(:,1:stimShift) = timeR(:,handleReg:handleReg+stimShift-1);
        elseif bhv.ResponseSide(iTrials) == 2
            rhandleChoiceR{iTrials}(:,1:stimShift) = timeR(:,handleReg:handleReg+stimShift-1);
        end
    end
    
    % get L/R choices as binary design matrix during stimulus and subsequent period.
    stimReg = round((preStimDur + stimTime(iTrials)) * sRate); %frame that is aligned with stimulus onset
    responseReg = ceil((preStimDur + spoutTime(iTrials)) * sRate);
    
    lstimChoiceR{iTrials} = false(size(timeR,1), frames);
    rstimChoiceR{iTrials} = false(size(timeR,1), frames);
    try
        temp = timeR(:,stimReg:responseReg-1);
        if bhv.ResponseSide(iTrials) == 1
            lstimChoiceR{iTrials}(:,1:responseReg-stimReg) = temp;
        elseif bhv.ResponseSide(iTrials) == 2
            rstimChoiceR{iTrials}(:,1:responseReg-stimReg) = temp;
        end
    end

    % get L/R choices as binary design matrix during response period. 2s max duration
    lresponseChoiceR{iTrials} = false(size(stimTemp,1),(2*sRate)+1);
    rresponseChoiceR{iTrials} = false(size(stimTemp,1),(2*sRate)+1);
    try
        if bhv.ResponseSide(iTrials) == 1
            lresponseChoiceR{iTrials}(:,1:frames-responseReg+1) = timeR(:,responseReg:end);
        elseif bhv.ResponseSide(iTrials) == 2
            rresponseChoiceR{iTrials}(:,1:frames-responseReg+1) = timeR(:,responseReg:end);
        end
    end
    
    % previous trial regressors
    if iTrials == 1 %don't use first trial
        prevRewardR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevChoiceR{iTrials} = NaN(size(timeR(:,1:end-4)));
        
    else %for all subsequent trials
        % same as for regular choice regressors but for prevoious trial
        prevChoiceR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.ResponseSide(bTrials(iTrials)-1) == 1
            prevChoiceR{iTrials} = timeR(:,1:end-4);
        end
        
        prevRewardR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.Rewarded(bTrials(iTrials)-1) %last trial was rewarded
            prevRewardR{iTrials} = timeR(:,1:end-4);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        rewardR{iTrials} = rewardR{iTrials}(:,1:gaussShift:end);
        prevRewardR{iTrials} = prevRewardR{iTrials}(:,1:gaussShift:end);
        lstimChoiceR{iTrials} = lstimChoiceR{iTrials}(:,1:gaussShift:end);
        rstimChoiceR{iTrials} = rstimChoiceR{iTrials}(:,1:gaussShift:end);
        lhandleChoiceR{iTrials} = lhandleChoiceR{iTrials}(:,1:gaussShift:end);
        rhandleChoiceR{iTrials} = rhandleChoiceR{iTrials}(:,1:gaussShift:end);
        lresponseChoiceR{iTrials} = lresponseChoiceR{iTrials}(:,1:gaussShift:end);
        rresponseChoiceR{iTrials} = rresponseChoiceR{iTrials}(:,1:gaussShift:end);
        prevChoiceR{iTrials} = prevChoiceR{iTrials}(:,1:gaussShift:end);
    end
    
    %determine timepoint of reward given
    waterR{iTrials} = false(frames, sRate * 2);
    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = round((preStimDur + water(iTrials)) * sRate); %timepoint in frames when reward was given
        if waterOn <= frames
            waterR{iTrials}(:, 1 : size(timeR,2) - waterOn + 1) = timeR(:, waterOn:end);
        end
    end
    
    if gaussShift > 1
        waterR{iTrials} = waterR{iTrials}(:,1:gaussShift:end); % subsample regressor
    end
    
    %% lever grabs
    cGrabs = levGrabL{iTrials};
    cGrabs(cGrabs >= postStimDur) = []; %remove grabs after end of imaging
    cGrabs(find(diff(cGrabs) < tapDur) + 1) = []; %remove grabs that are too close to one another
    lGrabR{iTrials} = histcounts(cGrabs,-preStimDur:1/sRate:postStimDur)'; %convert to binary trace
    
    cGrabs = levGrabR{iTrials};
    cGrabs(cGrabs >= postStimDur) = []; %remove grabs after end of imaging
    cGrabs(find(diff(cGrabs) < tapDur) + 1) = []; %remove grabs that are too close to one another
    rGrabR{iTrials} = histcounts(cGrabs,-preStimDur:1/sRate:postStimDur)'; %convert to binary trace
    
    %% pupil / whisk / nose / body regressors
    bhvFrameRate = round(1/mean(diff(pTime{bTrials(iTrials)}))); %framerate of face camera
    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab(iTrials) - preStimDur - pupilShift);
    trialTime = pTime{bTrials(iTrials)} - trialOn;
    rejIdx = trialTime < trialDur; %don't use late frames
    trialTime = trialTime(rejIdx);
    
    if isempty(trialTime) || trialTime(1) > 0 %check if there is missing time at the beginning of a trial
        warning(['Trial ' int2str(bTrials(iTrials)) ': Missing behavioral video frames at trial onset. Trial removed from analysis']);
        pupilR{iTrials} = NaN(frames, 1);
        whiskR{iTrials} = NaN(frames, 1);
        noseR{iTrials} = NaN(frames, 1);
        bodyR{iTrials} = NaN(frames, 1);
        
    else
        timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
        if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
            addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
            trialTime = [trialTime' addTime];
        end
        
        pupilR{iTrials} = Behavior_vidResamp(fPupil{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        pupilR{iTrials} = smooth(pupilR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        whiskR{iTrials} = Behavior_vidResamp(whisker{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        whiskR{iTrials} = smooth(whiskR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        noseR{iTrials} = Behavior_vidResamp(nose{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        noseR{iTrials} = smooth(noseR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        % body regressors
        bhvFrameRate = round(1/mean(diff(bTime{bTrials(iTrials)}))); %framerate of body camera
        trialTime = bTime{bTrials(iTrials)} - trialOn;
        rejIdx = trialTime < trialDur; %don't use late frames
        trialTime = trialTime(rejIdx);
        timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
        
        if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
            addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
            trialTime = [trialTime' addTime];
        end
        
        bodyR{iTrials} = Behavior_vidResamp(bodyM{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        bodyR{iTrials} = smooth(bodyR{iTrials}(end - frames + 1 : end), 'rlowess');
    end
    
    %% piezo sensor information
    if strcmpi(dType,'Widefield')
        if exist([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'file') ~= 2  %check if files exists on hdd and pull from server otherwise
            copyfile([sPath 'Analog_'  num2str(trials(iTrials)) '.dat'],[cPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
        end
        [~,Analog] = Widefield_LoadData([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'Analog'); %load analog data
        stimOn = find(diff(double(Analog(stimLine,:)) > 1500) == 1); %find stimulus onset in current trial
    elseif strcmpi(dType,'twoP')
        Analog = squeeze(data.analog(:,:,iTrials));
        stimOn = find(diff(double(Analog(stimLine,:)) > 1) == 1); %find stimulus onset in current trial
    end
    
    if ~isnan(stimTime(iTrials))
        try
            Analog(1,round(stimOn + ((postStimDur-stimTime(iTrials)) * 1000) - 1)) = 0; %make sure there are enough datapoints in analog signal
            temp = Analog(piezoLine,round(stimOn - ((preStimDur + stimTime(iTrials)) * 1000)) : round(stimOn + ((postStimDur - stimTime(iTrials))* 1000) - 1)); % data from piezo sensor. Should encode animals hindlimb motion.
            temp = smooth(double(temp), sRate*5, 'lowess')'; %do some smoothing
            temp = [repmat(temp(1),1,1000) temp repmat(temp(end),1,1000)]; %add some padding on both sides to avoid edge effects when resampling
            temp = resample(double(temp), sRate, 1000); %resample to imaging rate
            piezoR{iTrials} = temp(sRate + 1 : end - sRate)'; %remove padds again
            piezoR{iTrials} = piezoR{iTrials}(end - frames + 1:end); %make sure, the length is correct
        catch
            piezoR{iTrials} = NaN(frames, 1);
        end
    else
        piezoR{iTrials} = NaN(frames, 1);
    end
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
end

%% re-align behavioral video data and Vc to lever grab instead of stimulus onset
if strcmpi(dType,'Widefield')
    iiSpikeFrames = findInterictalSpikes(U, Vc, 2, false); %find interictal spikes
    Vc = interpOverInterictal(Vc, iiSpikeFrames); %interpolate over interictal spikes
end

V1 = reshape(V1, [], SessionData.nTrials, bhvDimCnt);
V2 = reshape(V2, [], SessionData.nTrials, bhvDimCnt);
dlcData = reshape(dlcData, [], SessionData.nTrials, length(dlcVars));

%if video sampling rate is different from widefield, resample video data
if bhvOpts.targRate ~= sRate
    vidR = NaN(size(Vc,2), length(bTrials), size(V1,3), 'single');
    moveR = NaN(size(Vc,2), length(bTrials), size(V1,3), 'single');
    dlcResamp = NaN(size(Vc,2), length(bTrials), size(dlcData,3), 'single');
    for iTrials = 1 : length(bTrials)
        temp = squeeze(V1(:,bTrials(iTrials),:));
        if ~any(isnan(temp(:)))
            trialTime = 1/bhvRate : 1/bhvRate : size(Vc,2)/sRate;
            vidR(:, iTrials, :) = Behavior_vidResamp(double(temp), trialTime, sRate);
            
            temp = squeeze(V2(:,bTrials(iTrials),:));
            moveR(:, iTrials, :) = Behavior_vidResamp(double(temp), trialTime, sRate);
            
            temp = squeeze(dlcData(:,bTrials(iTrials),:));
            dlcResamp(:, iTrials, :) = Behavior_vidResamp(double(temp), trialTime, sRate);
        end
    end
else
    %get correct trials from behavioral video data.
    vidR = V1(:,bTrials,:); 
    moveR = V2(:,bTrials,:);
    dlcResamp = dlcData(:,bTrials,:);
end
clear V1 V2 dlcData

% re-align video data
temp1 = NaN(dims,frames,trialCnt);
temp2 = NaN(frames,trialCnt,bhvDimCnt);
temp3 = NaN(frames,trialCnt,bhvDimCnt);
temp4 = NaN(frames,trialCnt,length(dlcVars));
temp5 = NaN(2,frames,trialCnt);
for x = 1 : size(vidR,2)
    try
        temp1(:,:,x) = Vc(:,(shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        temp2(:,x,:) = vidR((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        temp3(:,x,:) = moveR((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        temp4(:,x,:) = dlcResamp((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        if strcmpi(dType,'twoP')
            temp5(:,:,x) = data.DS(:,(shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        end
    catch
        fprintf(1,'Could not align trial %d. Relative stim time: %fs\n', x, stimTime(x));
    end
end
Vc = reshape(temp1,dims,[]); clear temp1
vidR = reshape(temp2,[],bhvDimCnt); clear temp2
moveR = reshape(temp3,[],bhvDimCnt); clear temp3
dlcResamp = reshape(temp4,[],length(dlcVars)); clear temp4

if strcmpi(dType,'twoP')
    DS = reshape(temp5,2,[]); %keep image motion trace for 2p imaging
end
clear temp5

%%  rebuild DLC data - filter and fill gaps in tracking
clear dlcLabels
dlcR = NaN(size(dlcResamp,1), size(dlcResamp,2) * 1/3, 'single');
Cnt = 0;
for x = 1 : 3 : length(dlcVars)
    
    temp = double(dlcResamp(:,x:x+1));
    rejIdx = dlcResamp(:,x+2) < 0.9; %dont use low-probability frames
    temp(rejIdx,:) = NaN;
    temp = fillgaps(temp, 100);

    dlcR(:, Cnt + (1:2)) = temp;
    dlcLabels(:, Cnt + (1:2)) = dlcVars(x:x+1);
    Cnt = Cnt + 2;
    
end
dlcR(isnan(dlcR)) = 0;
        
%% get proper design matrices for handle grab
lGrabR = cat(1,lGrabR{:});
lGrabR = Widefield_analogToDesign(lGrabR, 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift); %get design matrix
lGrabR = cat(1,lGrabR{:});

rGrabR = cat(1,rGrabR{:});
rGrabR = Widefield_analogToDesign(rGrabR, 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift); %get design matrix
rGrabR = cat(1,rGrabR{:});

% add DLC left paw labels
hIdx = contains(dlcLabels,'cam1_Elbow') | contains(dlcLabels,'cam1_Hand_L');
hTrace = zscore(dlcR(:,hIdx));
hTrace = smoothCol(hTrace, 1, 5, 'box');
[Q,~] = qr(zscore(hTrace),0); %orthogonalize
lGrabR = [lGrabR zscore(Q)]; %keep orthogonalized hand and elbow regressors

% add DLC right paw labels
hIdx = contains(dlcLabels,'cam1_Hand_R');
hTrace = zscore(dlcR(:,hIdx));
hTrace = smoothCol(hTrace, 1, 5, 'box');
[Q,~] = qr(zscore(hTrace),0); %orthogonalize
rGrabR = [rGrabR zscore(Q)]; %keep orthogonalized hand and elbow regressors

%% rebuild analog pupil regressor and get proper design matrices
temp = double(cat(1,pupilR{:}));
[b, a] = butter(4,0.1/sRate, 'high');
temp(~isnan(temp(:,1)),:) = single(filtfilt(b,a,temp(~isnan(temp(:,1)),:)));
temp = zscore(temp);

%find dilations
events = temp > pupilThresh;
events = [false(1,size(temp,2)); diff(events) == 1];
posMat = rateDisc_makeDesignMatrix(events, 3, opts);

%find constrictions
events = temp < -pupilThresh;
events = [false(1,size(temp,2)); diff(events) == 1];
negMat = rateDisc_makeDesignMatrix(events, ones(1,size(events,2))*3, opts);
pupilR = [single(temp),posMat{1},negMat{1}];

%% rebuild whisker regressors
[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,whiskR{:})), 1, trialCnt, sRate, sRate, motorIdx, gaussShift);
whiskR = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format

% get whiskers from cam1 and compute average
clear wTrace
wIdx = find(contains(dlcLabels,'cam1_Whisker'));
wTrace(:,1) = zscore(nanmean(dlcR(:, wIdx(1:2:length(wIdx))),2))'; %face whisker tracking
wTrace(:,2) = zscore(nanmean(dlcR(:, wIdx(2:2:length(wIdx))),2))'; %face whisker tracking
[b, a] = butter(4,0.1/sRate, 'high');
wTrace = single(filtfilt(b,a,double(wTrace)));

%keep an orthogonalized version of whisker position
traceOut = arrayfun(@(x) zscore(orthVec(whiskR(:,1), wTrace(:,x))), 1:size(wTrace,2), 'UniformOutput', false);
whiskR = [whiskR cat(2,traceOut{:})];

% compute differential and smoothed motion energy envelope
wTrace = smoothCol(wTrace, 1, 5, 'box');
wTrace = reshape(wTrace, [], trialCnt, 2);
wTrace = diff(wTrace,[],1);
wTrace = cat(1,wTrace(1,:,:), wTrace);
wTrace = reshape(wTrace, [] , 2);
wTrace = smoothCol(abs(hilbert(wTrace)), 1, 5, 'box');

% get event matrix
[dMat1, wTrace(:,1)] = Widefield_analogToDesign(double(wTrace(:,1)), dlcThresh, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, wTrace(:,2)] = Widefield_analogToDesign(double(wTrace(:,2)), dlcThresh, trialCnt, sRate, sRate, motorIdx, gaussShift);

%keep an orthogonalized version of absolute whisker velocity
traceOut = arrayfun(@(x) zscore(orthVec(whiskR(:,1), wTrace(:,x))), 1:size(wTrace,2), 'UniformOutput', false); 

% combine into larger matrix
whiskR = [whiskR, cat(2,traceOut{:}), cat(1,dMat1{:}), cat(1,dMat2{:})];


%% rebuild jaw regressors
clear jTrace
jIdx = find(contains(dlcLabels,'Jaw')); %this should give 4 regressors (2 from each camera)

%y-axis in cam1 and x-axis in cam2 should be very correlated. Average gives a good estimate of jaw opening/closing.
jTrace = dlcR(:, jIdx(2)).*-1; %invert cam1 tracking because of different camera orientation
jTrace = cat(2,jTrace,dlcR(:, jIdx(3)));
jTrace = zscore(nanmean(jTrace,2));
jawR = jTrace; %keep this trace for the model

% remove slow jaw movements by subtracting median in each trial
jTrace = reshape(double(jawR), [], trialCnt);
jTrace = bsxfun(@minus,jTrace,nanmedian(jTrace,1));
jTrace = jTrace(:);

% make design matrix
events = jTrace > dlcThresh;
events = [false(1,size(jTrace,2)); diff(events) == 1];
posMat = rateDisc_makeDesignMatrix(events, 3, opts);
jawR = [jawR posMat{1}]; %add design matrix

% add mouth regressors
clear mTrace
mIdx = find(contains(dlcLabels,'cam2_Mouth')); %this should give 4 regressors
%average x and y axis from left and right mout edge. They are very correlated so averaging will avoid degrading the model.
mTrace(:,1) = zscore(nanmean(dlcR(:, mIdx(1:2:length(mIdx))),2));
mTrace(:,2) = zscore(nanmean(dlcR(:, mIdx(2:2:length(mIdx))),2));

%orthogonalize against jaw movements. they are too correlated otherwise.
mTrace = arrayfun(@(x) zscore(orthVec(jawR(:,1), mTrace(:,x))), 1:size(mTrace,2), 'UniformOutput', false);

% add mouth regressors to jaw variable
jawR = [jawR cat(2,mTrace{:})];


%% rebuild nose regressors
noseR = cat(1,noseR{:});
noseR = zscore(smoothCol(noseR, 1, 5, 'box'));

% get whiskers from cam1 and compute average
clear nTrace
nIdx = contains(dlcLabels,'Nose') | contains(dlcLabels,'cam2_Whisker');
nTrace = zscore(dlcR(:, nIdx(2))); %nose tracking
nTrace = smoothCol(nTrace, 1, 5, 'box');
nTrace = reshape(nTrace, [], trialCnt);
nTrace = diff(nTrace,[],1);
nTrace = cat(1,nTrace(1,:), nTrace);
nTrace = nTrace(:);
nTrace = smoothCol(abs(hilbert(nTrace)), 1, 5, 'box');
[dMat, traceOut] = Widefield_analogToDesign(double(nTrace), dlcThresh, trialCnt, sRate, sRate, motorIdx, gaussShift);
noseR = [noseR zscore(traceOut) cat(1,dMat{:})];

% get other nose regressors
nTrace = zscore(dlcR(:, nIdx)); %nose tracking
[Q,~] = qr(cat(2,zscore(noseR(:,1:2)),nTrace),0); %orthogonalize
nTrace = zscore(Q(:,3:end));

% add DLC regressors to nose variable
noseR = [noseR nTrace];


%% rebuild body regressors
[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,bodyR{:})), 1, trialCnt, sRate, sRate, motorIdx, gaussShift);
bodyR = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format

shIdx = find(contains(dlcLabels,'cam2_Shoulder'));
shTrace = nanmean(dlcR(:,shIdx([1 3])),2); % add DLC shoulder regressors to body variable

% remove slow body movements by subtracting median in each trial
shTrace = reshape(double(shTrace), [], trialCnt);
shTrace = bsxfun(@minus,shTrace,nanmedian(shTrace,1));
shTrace = shTrace(:);
shTrace = (shTrace - nanmedian(shTrace)) ./ nanstd(shTrace);

% make design matrix
events = shTrace > dlcThresh;
events = [false(1,size(shTrace,2)); diff(events) == 1];
posMat = rateDisc_makeDesignMatrix(events, 3, opts);
bodyR = [bodyR shTrace posMat{1}]; %add design matrix

% include other shoulder DLC regressors
shTrace = zscore(dlcR(:,[2 4])); % add DLC shoulder regressors to body variable
[b, a] = butter(4,0.1/sRate, 'high');
shTrace = single(filtfilt(b,a,double(shTrace)));
bodyR = [bodyR shTrace]; %add design matrix


%% create ear variable
earTrace = zscore(nanmean(dlcR(:,contains(dlcLabels,'Ear')),2)); %average trace for ear position

% remove slow jaw movements by subtracting median in each trial
earTrace = reshape(double(earTrace), [], trialCnt);
earTrace = bsxfun(@minus,earTrace,nanmedian(earTrace,1));
earTrace = earTrace(:);

% get smoothed motion energy
earTrace = smoothCol(earTrace, 1, 5, 'box');
earTrace = reshape(earTrace, [], trialCnt);
earTrace = diff(earTrace,[],1);
earTrace = cat(1,earTrace(1,:), earTrace);
earTrace = earTrace(:);
earTrace = smoothCol(abs(hilbert(earTrace)), 1, 5, 'box');

% threshold and get design matrix
[dMat, earTrace] = Widefield_analogToDesign(earTrace, dlcThresh, trialCnt, sRate, sRate, motorIdx, gaussShift);
earR = [earTrace cat(1,dMat{:})];

% add dlc position signals
earTrace = zscore(dlcR(:,contains(dlcLabels,'Ear')));
[b, a] = butter(4,0.1/sRate, 'high');
earTrace = single(filtfilt(b,a,double(earTrace)));
earR = [earR single(earTrace)];
earR(isnan(earR)) = 0;

%% create grooming variable
groom = nanmean(dlcResamp(:,find(contains(dlcVars,'cam2'),1)+2:3: length(dlcVars)),2); %during grooming, marker probabilities drop down
events = groom < 0.9;
events = [false; diff(events) == 1];
negMat = rateDisc_makeDesignMatrix(events, 3, opts);
groomR = negMat{1};

%% create eye variable
eTrace = zscore(dlcR(:,contains(dlcLabels,'cam1_Eye'))); %this should give 4 regressors (2 from each camera)
eTrace = arrayfun(@(x) (removeOutlier(eTrace(:,x))),1:size(eTrace,2), 'UniformOutput', false);
eTrace = cat(2,eTrace{:});
eTrace = smoothCol(eTrace, 1, 5, 'box');
[Q,~] = qr(zscore(eTrace),0); %orthogonalize
eyeR = zscore(Q); %keep orthogonalized, analog eye regressors 

%% create licking variable
lLickR = cat(1,lLickR{:});
rLickR = cat(1,rLickR{:});

% get DLC labels to see premature licks
tIdx = find(contains(dlcVars,'cam2_Tongue')); %get tongue from bottom camera
tTrace = dlcResamp(:,tIdx(2)); % add DLC shoulder regressors to body variable

rejIdx = find(dlcResamp(:,tIdx(3)) < 0.9);
startAvg = mean(tTrace(rejIdx(diff(rejIdx) > 1) + 1)); %get all first confident measures to estimate the average
tTrace(rejIdx) = startAvg; %dont use low confidence markers
tTrace(isnan(tTrace)) = startAvg; %dont use low confidence markers
tTrace = smoothCol(tTrace, 1, 5, 'box');
tTrace = (tTrace - median(tTrace)) ./ std(tTrace);

%find right licks
events = tTrace > lickThresh;
events = [false(1,size(tTrace,2)); diff(events) == 1];
posMat = rateDisc_makeDesignMatrix(events, 3, opts);
lLickR = [lLickR posMat{1}];

%find left licks
events = tTrace < -lickThresh;
events = [false(1,size(tTrace,2)); diff(events) == 1];
negMat = rateDisc_makeDesignMatrix(events, ones(1,size(events,2))*3, opts);
rLickR = [rLickR negMat{1}];

% combine analog labels in lickR variable
lickR = [];
tIdx = find(contains(dlcVars,'Tongue')); %get tongue from bottom camera
for x = 1 : 3 : length(tIdx)
    rejIdx = find(dlcResamp(:,tIdx(3)) < 0.9);
    for y = 0:1
        tTrace = dlcResamp(:,tIdx(x+y)); % add DLC shoulder regressors to body variable
        startAvg = mean(tTrace(rejIdx(diff(rejIdx) > 1) + 1)); %get all first confident measures to estimate the average
        tTrace(rejIdx) = startAvg; %dont use low confidence markers
        tTrace(isnan(tTrace)) = startAvg; %dont use low confidence markers
        tTrace = smoothCol(tTrace, 1, 5, 'box');
        lickR = [lickR zscore(tTrace)];
    end
end
[Q,~] = qr(zscore(lickR),0); %orthogonalize
lickR = zscore(Q);
lickR(isnan(lickR)) = 0;

%% create hindpaw variable
[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,piezoR{:})), 1, trialCnt, sRate, sRate, motorIdx, gaussShift);
piezoR = [traceOut cat(1,dMat{:})]; %rebuild continuous format

%% reshape regressors, make design matrix and indices for regressors that are used for the model
leverInR = cat(1,leverInR{:});
leverInR(:,sum(leverInR) == 0) = [];

handleSoundR = cat(1,handleSoundR{:});

lTacStimR = cat(1,lTacStimR{:});
rTacStimR = cat(1,rTacStimR{:});
lAudStimR = cat(1,lAudStimR{:});
rAudStimR = cat(1,rAudStimR{:});

lfirstTacStimR = cat(1,lfirstTacStimR{:});
rfirstTacStimR = cat(1,rfirstTacStimR{:});
lfirstAudStimR = cat(1,lfirstAudStimR{:});
rfirstAudStimR = cat(1,rfirstAudStimR{:});

spoutR = cat(1,spoutR{:});
spoutOutR = cat(1,spoutOutR{:});
spoutR(:,sum(spoutR) == 0) = [];
spoutOutR(:,sum(spoutOutR) == 0) = [];

rewardR = cat(1,rewardR{:});
prevRewardR = cat(1,prevRewardR{:});

lstimChoiceR = cat(1,lstimChoiceR{:});
rstimChoiceR = cat(1,rstimChoiceR{:});
lhandleChoiceR = cat(1,lhandleChoiceR{:});
rhandleChoiceR = cat(1,rhandleChoiceR{:});
lresponseChoiceR = cat(1,lresponseChoiceR{:});
rresponseChoiceR = cat(1,rresponseChoiceR{:});
prevChoiceR = cat(1,prevChoiceR{:});
waterR = cat(1,waterR{:});

%% create full design matrix
fullR = [lhandleChoiceR rhandleChoiceR lstimChoiceR rstimChoiceR ...
    lresponseChoiceR rresponseChoiceR rewardR lGrabR rGrabR lLickR rLickR lickR ...
    handleSoundR lfirstTacStimR lTacStimR rfirstTacStimR rTacStimR ...
    lfirstAudStimR lAudStimR rfirstAudStimR rAudStimR prevRewardR prevChoiceR ...
    waterR piezoR whiskR noseR pupilR bodyR eyeR jawR groomR earR moveR vidR];

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels = [
    {'lhandleChoice'} {'rhandleChoice'} {'lstimChoice'} {'rstimChoice'} {'lresponseChoice'} {'rresponseChoice'}  ...
    {'reward'} {'lGrab'} {'rGrab'} {'lLick'} {'rLick'} {'lick'} {'handleSound'} ...
    {'lfirstTacStim'} {'lTacStim'} {'rfirstTacStim'} {'rTacStim'} {'lfirstAudStim'} {'lAudStim'} ...
    {'rfirstAudStim'} {'rAudStim'} {'prevReward'} {'prevChoice'} {'water'} {'piezo'} {'whisk'} ...
    {'nose'} {'pupil'} {'body'} {'eye'} {'jaw'} {'groom'} {'ear'} {'Move'} {'bhvVideo'}];

%index to reconstruct different response kernels
regIdx = [
    ones(1,size(lhandleChoiceR,2))*find(ismember(regLabels,'lhandleChoice')) ...
    ones(1,size(rhandleChoiceR,2))*find(ismember(regLabels,'rhandleChoice')) ...
    ones(1,size(lstimChoiceR,2))*find(ismember(regLabels,'lstimChoice')) ...
    ones(1,size(rstimChoiceR,2))*find(ismember(regLabels,'rstimChoice')) ...
    ones(1,size(lresponseChoiceR,2))*find(ismember(regLabels,'lresponseChoice')) ...
    ones(1,size(rresponseChoiceR,2))*find(ismember(regLabels,'rresponseChoice')) ...
    ones(1,size(rewardR,2))*find(ismember(regLabels,'reward')) ...
    ones(1,size(lGrabR,2))*find(ismember(regLabels,'lGrab')) ...
    ones(1,size(rGrabR,2))*find(ismember(regLabels,'rGrab')) ...
    ones(1,size(lLickR,2))*find(ismember(regLabels,'lLick')) ...
    ones(1,size(rLickR,2))*find(ismember(regLabels,'rLick')) ...
    ones(1,size(lickR,2))*find(ismember(regLabels,'lick')) ...
    ones(1,size(handleSoundR,2))*find(ismember(regLabels,'handleSound')) ...
    ones(1,size(lfirstTacStimR,2))*find(ismember(regLabels,'lfirstTacStim')) ...
    ones(1,size(lTacStimR,2))*find(ismember(regLabels,'lTacStim')) ...
    ones(1,size(rfirstTacStimR,2))*find(ismember(regLabels,'rfirstTacStim')) ...
    ones(1,size(rTacStimR,2))*find(ismember(regLabels,'rTacStim')) ...
    ones(1,size(lfirstAudStimR,2))*find(ismember(regLabels,'lfirstAudStim')) ...
    ones(1,size(lAudStimR,2))*find(ismember(regLabels,'lAudStim')) ...
    ones(1,size(rfirstAudStimR,2))*find(ismember(regLabels,'rfirstAudStim')) ...
    ones(1,size(rAudStimR,2))*find(ismember(regLabels,'rAudStim')) ...
    ones(1,size(prevRewardR,2))*find(ismember(regLabels,'prevReward')) ...
    ones(1,size(prevChoiceR,2))*find(ismember(regLabels,'prevChoice')) ...
    ones(1,size(waterR,2))*find(ismember(regLabels,'water')) ...
    ones(1,size(piezoR,2))*find(ismember(regLabels,'piezo')) ...
    ones(1,size(whiskR,2))*find(ismember(regLabels,'whisk')) ...
    ones(1,size(noseR,2))*find(ismember(regLabels,'nose')) ...
    ones(1,size(pupilR,2))*find(ismember(regLabels,'pupil')) ...
    ones(1,size(bodyR,2))*find(ismember(regLabels,'body')) ...
    ones(1,size(eyeR,2))*find(ismember(regLabels,'eye')) ...
    ones(1,size(jawR,2))*find(ismember(regLabels,'jaw')) ...
    ones(1,size(groomR,2))*find(ismember(regLabels,'groom')) ...
    ones(1,size(earR,2))*find(ismember(regLabels,'ear')) ...
    ones(1,size(moveR,2))*find(ismember(regLabels,'Move')) ...
    ones(1,size(vidR,2))*find(ismember(regLabels,'bhvVideo'))];

% orthogonalize video against spout/handle movement
vidIdx = ismember(regIdx, find(ismember(regLabels,{'Move' 'bhvVideo'}))); %index for video regressors
trialIdx = ~isnan(mean(fullR(:,vidIdx),2)); %don't use trials that failed to contain behavioral video data

smallR = [leverInR spoutR spoutOutR lGrabR rGrabR lLickR rLickR lickR piezoR whiskR noseR ...
          pupilR bodyR eyeR jawR groomR earR moveR vidR];

[Q, ~] = qr(smallR(trialIdx,:),0);
vidR = fullR(:,vidIdx); %keep original video regressors
fullR(trialIdx,vidIdx) = zscore(Q(:, end - sum(vidIdx) + 1 : end));

% reject trials with broken regressors that contain NaNs
trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fprintf(1, 'Rejected %d/%d trials for NaN entries in regressors\n', sum(trialIdx)/frames, trialCnt);
fullR(trialIdx,:) = []; %clear bad trials
vidR(trialIdx,:) = []; %clear bad trials

%% run QR and check for rank-defficiency
rejIdx = nansum(abs(fullR)) < 10;
[~, fullQRR] = qr(bsxfun(@rdivide,fullR(:,~rejIdx),sqrt(sum(fullR(:,~rejIdx).^2))),0); %orthogonalize design matrix
if ispc
figure; plot(abs(diag(fullQRR))); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
end
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < sum(~rejIdx) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end

% reject regressors that are too sparse or rank-defficient
fullR(:,rejIdx) = []; %clear empty regressors
fprintf(1, 'Rejected %d/%d empty regressors\n', sum(rejIdx),length(rejIdx));

%% save modified Vc
Vc(:,trialIdx) = []; %clear bad trials
Vc = bsxfun(@minus, Vc, mean(Vc,2)); %should be zero-mean

if strcmpi(dType,'Widefield')
    save([cPath 'interpVc.mat'], 'Vc', 'frames', 'preStimDur', 'postStimDur', 'bTrials');
elseif strcmpi(dType,'twoP')
    DS(:,trialIdx) = []; %clear bad trials
    save([cPath 'interpVc.mat'], 'Vc', 'DS', 'frames', 'preStimDur', 'postStimDur');
end

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
clear stimR lGrabR rGrabR waterR lLickR rLickR lAudStimR rAudStimR rewardR prevRewardR ...
    lstimChoiceR rstimChoiceR lhandleChoiceR rhandleChoiceR ...
    prevChoiceR pupilR moveR piezoR whiskR noseR bodyR

%% run ridge regression in low-D
% run model. Zero-mean without intercept. only video qr.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'orgdimBeta.mat'], 'dimBeta', 'ridgeVals');
save([cPath filesep 'orgregData.mat'], 'fullR', 'spoutR', 'leverInR', 'rejIdx' ,'trialIdx', 'regIdx', 'regLabels','gaussShift','fullQRR','-v7.3');
[Vm, fullBeta, fullR, fullIdx, fullRidge, fullLabels, fullMap, fullMovie] = crossValModel(regLabels);
save([cPath 'orgfullcorr.mat'], 'Vm', 'fullBeta', 'fullIdx', 'fullR', 'fullLabels', 'fullRidge', 'regLabels', 'fullMap', 'fullMovie','-v7.3');

mInd = ismember(regIdx(~rejIdx), find(ismember(regLabels,motorLabels)));
motorR = fullR(:, mInd);
[motorRidge, motorBeta] = ridgeMML(Vc', motorR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, motor only model: %f\n', mean(motorRidge));
Vm = (motorR * motorBeta)';
save([cPath 'orgVMotor.mat'], 'Vm', 'frames'); %save predicted data based on motor model

mInd = ismember(regIdx(~rejIdx), find(ismember(regLabels,motorLabels(~ismember(motorLabels,opMotorLabels)))));
spontMotorR = fullR(:, mInd);
[spontMotorRidge, spontMotorBeta] = ridgeMML(Vc', spontMotorR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, spont-motor only model: %f\n', mean(spontMotorRidge));
Vm = (spontMotorR * spontMotorBeta)';
save([cPath 'orgVspontMotor.mat'], 'Vm', 'frames'); %save predicted data based on spont motor model

% run video model
[vidRidge, vidBeta] = ridgeMML(Vc', vidR, true); %get ridge penalties and beta weights.
save([cPath 'orgVidBeta.mat'], 'vidBeta', 'vidRidge');
save([cPath filesep 'vidRegData.mat'], 'vidR', 'trialIdx', 'frames', '-v7.3');

%% run motor/task/opMotor and spontMotor only models. Zero-mean without intercept.
cIdx = ismember(regIdx(~rejIdx), find(ismember(regLabels,motorLabels))); %get index for motor regressors
motorLabels = regLabels(sort(find(ismember(regLabels,motorLabels)))); %make sure motorLabels is in the right order

[Vmotor, motorBeta, motorR, motorIdx, motorRidge, motorLabels, motorMap, motorMovie] = crossValModel(motorLabels);
fprintf('Mean ridge penalty for motor-only, zero-mean model: %f\n', mean(motorRidge));
save([cPath 'interpVmotor.mat'], 'Vmotor', 'frames'); %save predicted data based on motor model
save([cPath 'motorBeta.mat'], 'motorBeta', 'motorRidge');
save([cPath filesep 'motorregData.mat'], 'motorR','trialIdx', 'motorIdx', 'motorLabels', 'motorMap', 'motorMovie', 'gaussShift','-v7.3');

[Vtask, taskBeta, taskR, taskIdx, taskRidge, taskLabels, taskMap, taskMovie] = crossValModel(regLabels(~ismember(regLabels,motorLabels)));
fprintf('Mean ridge penalty for task-only, zero-mean model: %f\n', mean(taskRidge));
save([cPath 'interpVtask.mat'], 'Vtask', 'frames'); %save predicted data based on task model
save([cPath 'taskBeta.mat'], 'taskBeta', 'taskRidge');
save([cPath filesep 'taskregData.mat'], 'taskR','trialIdx', 'taskIdx', 'taskLabels', 'taskMap', 'taskMovie', 'gaussShift','-v7.3');

[VtaskOpMotor, taskOpMotorBeta, taskOpMotorR, taskOpMotorIdx, taskOpMotorRidge, taskOpMotorLabels, taskOpMotorMap, taskOpMotorMovie] = crossValModel([regLabels(~ismember(regLabels,motorLabels)) opMotorLabels]);
fprintf('Mean ridge penalty for task+opMotor, zero-mean model: %f\n', mean(taskOpMotorRidge));
save([cPath 'interpVtaskOpMotor.mat'], 'VtaskOpMotor', 'frames'); %save predicted data based on taskOpMotor model
save([cPath 'taskOpMotorBeta.mat'], 'taskOpMotorBeta', 'taskOpMotorRidge');
save([cPath filesep 'taskOpMotorregData.mat'], 'taskOpMotorR','trialIdx', 'taskOpMotorIdx', 'taskOpMotorLabels', 'taskOpMotorMap', 'taskOpMotorMovie', 'gaussShift','-v7.3');

[VtaskSpMotor, taskSpMotorBeta, taskSpMotorR, taskSpMotorIdx, taskSpMotorRidge, taskSpMotorLabels, taskSpMotorMap, taskSpMotorMovie] = crossValModel(regLabels(~ismember(regLabels,opMotorLabels)));
fprintf('Mean ridge penalty for task+opMotor, zero-mean model: %f\n', mean(taskSpMotorRidge));
save([cPath 'interpVtaskSpMotor.mat'], 'VtaskSpMotor', 'frames'); %save predicted data based on taskSpMotor model
save([cPath 'taskSpMotorBeta.mat'], 'taskSpMotorBeta', 'taskSpMotorRidge');
save([cPath filesep 'taskSpMotorregData.mat'], 'taskSpMotorR','trialIdx', 'taskSpMotorIdx', 'taskSpMotorLabels', 'taskSpMotorMap', 'taskSpMotorMovie', 'gaussShift','-v7.3');

[VopMotor, opMotorBeta, opMotorR, opMotorIdx, opMotorRidge, opMotorLabels, opMotorMap, opMotorMovie] = crossValModel(opMotorLabels);
fprintf('Mean ridge penalty for opMotor-only, zero-mean model: %f\n', mean(opMotorRidge));
save([cPath 'interpVopMotor.mat'], 'VopMotor', 'frames'); %save predicted data based on operant motor model
save([cPath 'opMotorBeta.mat'], 'opMotorBeta', 'opMotorRidge');
save([cPath filesep 'opMotorregData.mat'], 'opMotorR','trialIdx', 'opMotorIdx', 'opMotorLabels', 'opMotorMap', 'opMotorMovie', 'gaussShift','-v7.3');

[VspontMotor, spontMotorBeta, spontMotorR, spontMotorIdx, spontMotorRidge, spontMotorLabels, spontMotorMap, spontMotorMovie] = crossValModel(motorLabels(~ismember(motorLabels,opMotorLabels)));
fprintf('Mean ridge penalty for spontMotor-only, zero-mean model: %f\n', mean(spontMotorRidge));
save([cPath 'interpVspontMotor.mat'], 'VspontMotor', 'frames'); %save predicted data based on spontaneous motor model
save([cPath 'spontMotorBeta.mat'], 'spontMotorBeta', 'spontMotorRidge');
save([cPath filesep 'spontMotorregData.mat'], 'spontMotorR','trialIdx', 'spontMotorIdx', 'spontMotorLabels', 'spontMotorMap', 'spontMotorMovie', 'gaussShift','-v7.3');


%% run DLC model for comparison
% rebuild dlcR again
clear dlcLabels
dlcR = NaN(size(dlcResamp,1), size(dlcResamp,2) * 1/3);
Cnt = 0;
for x = 1 : 3 : length(dlcVars)
    temp = double(dlcResamp(:,x:x+1));
    temp(dlcResamp(:,x+2) < 0.9,:) = NaN; %dont use low-probability frames
    temp = fillgaps(temp, 100);
    
    dlcR(:, Cnt + (1:2)) = temp;
    dlcLabels(:, Cnt + (1:2)) = dlcVars(x:x+1);
    Cnt = Cnt + 2;
end

% filter and smooth
[b, a] = butter(4,0.1/sRate, 'high');
dlcR = single(filtfilt(b,a,dlcR));
dlcR = smoothCol(single(dlcR), 1, 5, 'box');
dlcR(trialIdx, :) = [];

% check DLC model
[dlcRidge, dlcBeta] = ridgeMML(Vc', dlcR, true); %get ridge penalties and beta weights.
save([cPath 'dlcBeta.mat'], 'dlcBeta', 'dlcRidge');

regIdx = [regIdx ones(1,size(dlcR,2)).*(max(regIdx)+1)];
rejIdx = [rejIdx, false(1,size(dlcR,2))];
regLabels = [regLabels dlcLabels]; fullR = [fullR dlcR];
[Vdlc, dlcBeta, dlcR, dlcIdx, dlcRidge, dlcLabels, dlcMap, dlcMovie] = crossValModel(dlcLabels);
save([cPath filesep 'dlcregData.mat'], 'dlcR','trialIdx', 'dlcLabels', 'gaussShift','-v7.3');


% %% subtract mean from Vc and run model again
% Vc = reshape(Vc, size(Vc,1), frames, []);
% meanV = mean(Vc,3);
% for iTrials = 1 : size(Vc,3)
%     Vc(:,:,iTrials) = bsxfun(@minus,Vc(:,:,iTrials),meanV); %subtract trial-average
% end
% Vc = reshape(Vc, size(Vc,1), []);
% Vc = bsxfun(@minus, Vc, mean(Vc,2)); %should be zero-mean
% 
% if strcmpi(dType,'Widefield')
%     save([cPath 'nomeanVc.mat'], 'Vc', 'meanV', 'frames', 'preStimDur', 'postStimDur');
% elseif strcmpi(dType,'twoP')
%     save([cPath 'nomeanVc.mat'], 'Vc', 'DS', 'frames', 'preStimDur', 'postStimDur');
% end
% 
% % re-run model. Zero-mean without intercept.
% [ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
% fprintf('Mean ridge penalty for mean-subtracted Vc with zero-mean model: %f\n', mean(ridgeVals));
% save([cPath 'nomeandimBeta.mat'], 'dimBeta', 'ridgeVals');
% save([cPath filesep 'nomeanregData.mat'], 'fullR', 'spoutR', 'leverInR', 'rejIdx' ,'trialIdx', 'regIdx', 'regLabels','gaussShift','fullQRR','-v7.3');
% [~, ~, ~, ~, ~, fullLabels, fullMap, fullMovie] = crossValModel(regLabels);
% save([cPath 'nomeanfullcorr.mat'], 'fullLabels', 'fullMap', 'fullMovie','-v7.3');
% 
% [~, ~, ~, ~, motorRidge, motorLabels, motorMap, motorMovie] = crossValModel(motorLabels);
% fprintf('Mean ridge penalty for motor-only, zero-mean model: %f\n', mean(motorRidge));
% save([cPath filesep 'nomeanmotorregData.mat'], 'motorLabels', 'motorMap', 'motorMovie','-v7.3');
% 
% [~, ~, ~, ~, taskRidge, taskLabels, taskMap, taskMovie] = crossValModel(regLabels(~ismember(regLabels,motorLabels)));
% fprintf('Mean ridge penalty for task-only, zero-mean model: %f\n', mean(taskRidge));
% save([cPath filesep 'nomeantaskregData.mat'], 'taskLabels', 'taskMap', 'taskMovie', '-v7.3');



%% nested functions
    function [Vm, cBeta, cR, subIdx, cRidge, cLabels, cMap, cMovie] =  crossValModel(cLabels)
        
        cIdx = ismember(regIdx(~rejIdx), find(ismember(regLabels,cLabels))); %get index for task regressors
        cLabels = regLabels(sort(find(ismember(regLabels,cLabels)))); %make sure motorLabels is in the right order
        
        %create new regressor index that matches motor labels
        subIdx = regIdx(~rejIdx);
        subIdx = subIdx(cIdx);
        temp = unique(subIdx);
        for x = 1 : length(temp)
            subIdx(subIdx == temp(x)) = x;
        end
        cR = fullR(:,cIdx);
        
        Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
        randIdx = randperm(size(Vc,2)); %generate randum number index
        foldCnt = floor(size(Vc,2) / ridgeFolds);
        cBeta = cell(1,ridgeFolds);
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(Vc,2));
            
            if ridgeFolds > 1
                dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
                if iFolds == 1
                    [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
                else
                    [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
                end
                Vm(:,~dataIdx) = (cR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data
                
                if rem(iFolds,ridgeFolds/5) == 0
                    fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
                end
            else
                [cRidge, cBeta{iFolds}] = ridgeMML(Vc', cR, true); %get beta weights for task-only model.
                Vm = (cR * cBeta{iFolds})'; %predict remaining data
                disp('Ridgefold is <= 1, fit to complete dataset instead');
            end
        end
        
        % computed all predicted variance
        Vc = reshape(Vc,size(Vc,1),[]);
        Vm = reshape(Vm,size(Vm,1),[]);
        if length(size(U)) == 3
            U = arrayShrink(U, squeeze(isnan(U(:,:,1))));
        end
        covVc = cov(Vc');  % S x S
        covVm = cov(Vm');  % S x S
        cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
        covP = sum((U * cCovV) .* U, 2)';  % 1 x P
        varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
        varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
        stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
        cMap = gather((covP ./ stdPxPy)');
        
        % movie for predicted variance
        cMovie = zeros(size(U,1),frames, 'single');
        for iFrames = 1:frames
            
            frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
            cData = bsxfun(@minus, Vc(:,frameIdx), mean(Vc(:,frameIdx),2));
            cModel = bsxfun(@minus, Vm(:,frameIdx), mean(Vm(:,frameIdx),2));
            covVc = cov(cData');  % S x S
            covVm = cov(cModel');  % S x S
            cCovV = cModel * cData' / (length(frameIdx) - 1);  % S x S
            covP = sum((U * cCovV) .* U, 2)';  % 1 x P
            varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
            varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
            stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
            cMovie(:,iFrames) = gather(covP ./ stdPxPy)';
            clear cData cModel
            
        end
        fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));
        
        
    end
end

