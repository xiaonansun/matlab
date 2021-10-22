function newVc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts, varargin)
% code to re-align Vc so each trial is aligned to different task episodes.
% Alignment is done to baseline, handle, stimulus and delay period.
% 'segFrames' defines how many frames per task episodes should be in the
% newVc output. SegFrames should be cumulative, so say 'baseline should be
% from trial 1 to segFrames(1), handle data should be from segFrames(1) to
% segFrames(2) and so on...
% 2021-09-26 modified by Richard: two optional inputs can specify the
% animal and session ID, which in turn allows Vc to be saved
% 2021-10-04 added some default parameters so the function can run without
% segFrames, opts, animal, and session

trialStimFrame = 93; 
msPerFrame = 32.3638;
segIdx = [1 0.75 1.25 0.5 1];

if ~exist('opts','var') || isempty(opts)
    opts.preStim = trialStimFrame*msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
    opts.frameRate = 1000/msPerFrame; % Frame rate of imaging
    sRate = opts.frameRate;
end

if ~exist('segFrames','var') || isempty(segFrames)
    segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment
end

if nargin == 0
    disp('Missing animal and session ID imnput, data will not be saved.')
elseif nargin == 1
    disp('Missing session input, data will not be saved.')
elseif nargin == 2
    animal = varargin{1}; session = varargin{2};
end

S = twoP_settings;

%% align imaging data using Session data
rejCnt = 0;
newVc = NaN(size(Vc,1), segFrames(5), size(Vc,3), 'single'); %new Vc to capture max duration of each segment
for iTrials = 1 : size(Vc,3)
    
    % get indices for current trial
    stimOn = cBhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    handleOn = [reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    clear cIdx
    cIdx(1) = handleOn(find(handleOn == cBhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    cIdx(2) = stimOn;
    cIdx(3) = max(cat(2,cBhv.stimEvents{iTrials}{:})) + stimOn; %time of last stimulus event
    cIdx(4) = cBhv.RawEvents.Trial{iTrials}.States.MoveSpout(1);
    cIdx = floor((cIdx - stimOn + opts.preStim) * opts.frameRate); %convert to frames. This is the last frame of each segment.
    cIdx(end + 1) = size(Vc,2);
    
    if cIdx(1) > 0 %in very rare cases there might be something wrong with handle time. don't use those trials.
        if segFrames(1) >= cIdx(1)
            newVc(:, segFrames(1) - cIdx(1) + 1 : segFrames(1), iTrials) = Vc(:, 1 : cIdx(1), iTrials); % baseline
        elseif segFrames(1) < cIdx(1)
            newVc(:, 1 : segFrames(1), iTrials) = Vc(:, cIdx(1) - segFrames(1) + 1 : cIdx(1), iTrials); % baseline
        end
           
        newVc(:, segFrames(1) + 1 : segFrames(1) + (diff(cIdx(1:2))), iTrials) = Vc(:, cIdx(1) + 1 : cIdx(2), iTrials); %handle period
        newVc(:, segFrames(2) + 1 : segFrames(2) + (diff(cIdx(2:3))), iTrials) = Vc(:, cIdx(2) + 1 : cIdx(3), iTrials); %stimulus period
        
        maxDiff = min([segFrames(3) + (diff(cIdx(3:4))) segFrames(5)]) - segFrames(3); %maximal possible delay duration
        newVc(:, segFrames(3) + 1 : segFrames(3) + maxDiff, iTrials) = Vc(:, cIdx(3) + 1 : cIdx(3) + maxDiff, iTrials); %delay period
        
        if segFrames(4) + (diff(cIdx(4:5))) > segFrames(5)
            newVc(:, segFrames(4) + 1 : segFrames(5), iTrials) = Vc(:, cIdx(4) + 1 : cIdx(4) + (segFrames(5) - segFrames(4)), iTrials); %response period
        else
            newVc(:, segFrames(4) + 1 : segFrames(4) + (diff(cIdx(4:5))), iTrials) = Vc(:, cIdx(4) + 1 : cIdx(5), iTrials); %response period
        end
    else
        rejCnt = rejCnt + 1;
    end
end

if rejCnt > 0
    warning(['!!! Couldnt use ' num2str(rejCnt) ' trials because of broken handle initialization time !!!'])
end

%% Save Vc
Vc = newVc;

if exist('animal','var') && exist('session','var')
    VcSavePath = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'Vc.mat');
    save(VcSavePath,'Vc');
    disp(['Saved Vc as ' VcSavePath]);
end