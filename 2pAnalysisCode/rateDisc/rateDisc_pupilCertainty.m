function [bhv, allPupil, allTimes] = rateDisc_pupilCertainty(fPath,baseDur)
% code to check if mouse pupil shows features of decision uncertainty

% baseDur = 3;          % Duration of baseline before lever grab in seconds
trialDur = 8;         % Duration of a trial in seconds
sRate = 30;           % video framerate

%%
bhvFile = dir([fPath '*SpatialDisc*.mat']);
load([fPath filesep bhvFile.name],'SessionData')
load([fPath filesep 'Vc.mat'],'bTrials')
load([fPath filesep 'BehaviorVideo' filesep 'FilteredPupil.mat'], 'fillPupil', 'pTime')

SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dat

timeCheck1 = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
if (timeCheck1 > 3590 && timeCheck1 < 3610) %timeshift by one hour (+- 10seconds)
    warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
    for iTrials = 1 : length(pTime)
        pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
    end
elseif timeCheck1 > 30 || timeCheck1 < -30
    error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
end

%%
stimTime = NaN(1,length(bTrials));
stimEnd = NaN(1,length(bTrials));
spoutTime = NaN(1,length(bTrials));
allPupil = NaN(trialDur*sRate, length(bTrials));

for iTrials = 1 : length(bTrials)
    
    cPupil = smooth(fillPupil{bTrials(iTrials)}, 15, 'rlowess');
    
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    stimGrab = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset

    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab - baseDur);
    trialTime = pTime{bTrials(iTrials)} - trialOn;
    firstFrame = find(trialTime > 0, 1);
    
    stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab + baseDur; %time of stimulus onset - measured from soundcard
    stimTime(iTrials) = round(stimTime(iTrials) * sRate);
    stimEnd(iTrials) = stimTime(iTrials) + round(max([bhv.stimEvents{iTrials}{:}]) * sRate);
    
    spoutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - stimGrab + baseDur;
    spoutTime(iTrials) = round(spoutTime(iTrials) * sRate);
    
    try
        cPupil = cPupil(firstFrame : firstFrame + (trialDur*sRate) - 1);
    catch
        cPupil = cPupil(firstFrame:end);
    end
    allPupil(1:length(cPupil(firstFrame:end)),iTrials) = cPupil(firstFrame:end);
end

allTimes = [stimTime; stimEnd; spoutTime];
    
%%

% 
% function [bhv, allPupil, evidenceStrength] = rateDisc_pupilCertainty(fPath)
% % code to check if mouse pupil shows features of decision uncertainty
% 
% baseDur = 3;          % Duration of baseline before lever grab in seconds
% baseWin = 0.5;        % time used for averaging
% stimWin = 0.5;        % time used for averaging
% delayWin = 0.5;       % time used for averaging
% respWin = 0.5;        % time used for averaging
% pupilShift = 0.5;     % time shift of pupil trace to account for muscle
% 
% %%
% bhvFile = dir([fPath '*SpatialDisc*.mat']);
% load([fPath filesep bhvFile.name],'SessionData')
% load([fPath filesep 'Vc.mat'],'bTrials')
% load([fPath filesep 'BehaviorVideo' filesep 'FilteredPupil.mat'], 'fillPupil', 'pTime')
% 
% SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds
% bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dat
% 
% timeCheck1 = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
% if (timeCheck1 > 3590 && timeCheck1 < 3610) && (timeCheck2 > 3590 && timeCheck2 < 3610) %timeshift by one hour (+- 10seconds)
%     warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
%     for iTrials = 1 : length(pTime)
%         pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
%     end
% elseif timeCheck1 > 30 || timeCheck1 < -30
%     error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
% end
% 
% %%
% allPupil = NaN(baseWin, length(bTrials));
% evidenceStrength = NaN(length(bTrials),1);
% 
% for iTrials = 1 : length(bTrials)
%     
%     cPupil = smooth(fillPupil{bTrials(iTrials)}, 15, 'rlowess');
%     
%     leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
%         reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
%         reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
%     
%     stimGrab = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
%     bhvFrameRate = round(1/mean(diff(pTime{bTrials(iTrials)}))); %framerate of face camera
%     trialOn = bhv.TrialStartTime(iTrials) + (stimGrab - baseDur);
%     trialTime = pTime{bTrials(iTrials)} - trialOn;
%     
%     stimTime = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab + baseDur + pupilShift; %time of stimulus onset - measured from soundcard
%     stimEnd = stimTime + max([bhv.stimEvents{iTrials}{:}]); %time of last stimulus event
%     spoutTime = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - stimGrab + baseDur + pupilShift;
%     
%     evidenceStrength(iTrials) = abs(sum(bhv.stimEvents{iTrials}{1} < 1) - sum(bhv.stimEvents{iTrials}{2} < 1));
%     
%     %
%     baseIdx = trialTime > 0 & trialTime < baseWin;
%     meanP(iTrials,1) = mean(cPupil(baseIdx));
%     dataLength(iTrials,1) = sum(baseIdx);
%     
%     stimIdx = trialTime > stimTime & trialTime < stimTime + stimWin;
%     stimIdx(trialTime > stimEnd) = false; %dont use times when stimulus is already over
%     meanP(iTrials,2) = mean(cPupil(stimIdx));
%     dataLength(iTrials,2) = sum(stimIdx);
%     
%     delayIdx = trialTime > stimEnd & trialTime < stimEnd + delayWin;
%     delayIdx(trialTime > spoutTime) = false; %dont use times when stimulus is already over
%     meanP(iTrials,3) = mean(cPupil(delayIdx));
%     dataLength(iTrials,3) = sum(delayIdx);
%     
%     respIdx = trialTime > spoutTime & trialTime < spoutTime + respWin;
%     meanP(iTrials,4) = mean(cPupil(respIdx));
%     dataLength(iTrials,4) = sum(respIdx);
%     
%     %
%     firstFrame = find(trialTime > 0, 1);
%     allPupil(1:length(cPupil(firstFrame:end)),iTrials) = cPupil(firstFrame:end);
%     
% end
