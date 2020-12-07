fPath = 'Z:\data\Behavior_Simon\CSP32\SpatialDisc\Session Data\CSP32_SpatialDisc_Jul23_2020_Session2.mat';
load(fPath);
bhv=SessionData;
trialCnt=bhv.cTrial;

%%
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

%% Compute the duration between stimulus onset and the second lick (committed lick)

minLicks=2;
a=find(cellfun(@numel,lickL)>=2); b=find(cellfun(@numel,lickR)>=2);
aa=cellfun(@diff,lickL,'UniformOutput',false); aaa=cell2mat(aa); 
bb=cellfun(@diff,lickR,'UniformOutput',false); bbb=cell2mat(bb);

secondLick

figure(1); histogram(aaa,100); hold on;
histogram(bbb,100);
% maxSpoutRegs = length(min(round((preStimDur + spoutTime) * sRate)) : frames); %maximal number of required spout regressors