function newVc = rateDisc_getRealignment(fullR, Vc, regIdx, regLabels, segFrames, frames)
%% reshape data into trials and get indices
cSegFrames = cumsum(segFrames);
Vc = reshape(Vc, size(Vc,1), frames, []); 
fullR = reshape(fullR, frames, [], size(fullR,2));
newVc = NaN(size(Vc,1), sum(segFrames), size(Vc,3), 'single'); %new Vc to capture max duration of each segment

segIdx = [find(regIdx == find(ismember(regLabels,'lhandleChoice')), 1) ...
    find(regIdx == find(ismember(regLabels,'rhandleChoice')), 1) ...
    find(regIdx == find(ismember(regLabels,'lstimChoice')), 1) ...
    find(regIdx == find(ismember(regLabels,'rstimChoice')), 1) ...
    find(regIdx == find(ismember(regLabels,'lresponseChoice')), 1) ...
    find(regIdx == find(ismember(regLabels,'rresponseChoice')), 1)];

stimIdx = [find(regIdx == find(ismember(regLabels,'lTacStim')), 1) ...
    find(regIdx == find(ismember(regLabels,'rTacStim')), 1) ...
    find(regIdx == find(ismember(regLabels,'lAudStim')), 1) ...
    find(regIdx == find(ismember(regLabels,'rAudStim')), 1)];

%% baseline should be constant across trials
handleOn = find(nanmean(sum(fullR(:,:,segIdx(1:2)),3),2));
if cSegFrames(1) > handleOn
    newVc(:, cSegFrames(1) - handleOn + 1 : cSegFrames(1), :) = Vc(:, 1 : handleOn, :); %transfer baseline    
elseif cSegFrames(1) <= handleOn
    newVc(:, 1 : cSegFrames(1), :) = Vc(:, handleOn - cSegFrames(1) + 1 : handleOn, :); %transfer baseline
end
           
%% loop through trials and align subsequent periods
for iTrials = 1 : size(Vc,3)
    
    stimOn = find(sum(fullR(:,iTrials,segIdx(3:4)),3)); %index for stimulus onset
    newVc(:, cSegFrames(1) + 1 : stimOn, iTrials) = Vc(:, cSegFrames(1) + 1 : stimOn, iTrials); %handle period
    
    stimOff = find(sum(fullR(:,iTrials, stimIdx), 3), 1, 'last'); % last frame of stimulus sequence
    newVc(:, cSegFrames(2) + 1 : cSegFrames(2) + (stimOff - stimOn), iTrials) = Vc(:, stimOn + 1 : stimOff, iTrials); %stimulus period

    spoutIn = find(sum(fullR(:,iTrials,segIdx(5:6)),3)); %end of delay period
    newVc(:, cSegFrames(3) + 1 : cSegFrames(3) + (spoutIn - stimOff), iTrials) = Vc(:, stimOff + 1 : spoutIn, iTrials); %delay period

    if cSegFrames(4) + (size(Vc,2) - spoutIn) > cSegFrames(5)
        newVc(:, cSegFrames(4) + 1 : cSegFrames(5), iTrials) = Vc(:, spoutIn + 1 : spoutIn + diff(cSegFrames(4:5)), iTrials); %response period
    else
        newVc(:, cSegFrames(4) + 1 : cSegFrames(4) + (size(Vc,2) - spoutIn), iTrials) = Vc(:, spoutIn + 1 : end, iTrials); %response period
    end
end

%% alternate approach using Session data
% newVc = NaN(size(Vc,1), segFrames(5), size(Vc,3), 'single'); %new Vc to capture max duration of each segment
% bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset
% 
% for iTrials = 1 : size(Vc,3)
% 
%     % get indices for current trial
%     stimOn = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
%     handleOn = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
%         reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
%         reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
%     
%     clear cIdx
%     cIdx(1) = handleOn(find(handleOn == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
%     cIdx(2) = stimOn;
%     cIdx(3) = max(cat(2,bhv.stimEvents{iTrials}{:})) + stimOn; %time of last stimulus event
%     cIdx(4) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1);
%     cIdx = floor((cIdx - stimOn + opts.preStim) * opts.frameRate); %convert to frames. This is the last frame of each segment.
%     cIdx(end + 1) = size(Vc,2);
%     
%     if segFrames(1) > cIdx(1)
%         newVc(:, segFrames(1) - cIdx(1) + 1 : segFrames(1), iTrials) = Vc(:, 1 : cIdx(1), iTrials); % baseline
%     elseif segFrames(1) < cIdx(1)
%         newVc(:, 1 : segFrames(1), iTrials) = Vc(:, cIdx(1) - segFrames(1) + 1 : cIdx(1), iTrials); % baseline
%     end
%         
%     newVc(:, segFrames(1) + 1 : segFrames(1) + (diff(cIdx(1:2))), iTrials) = Vc(:, cIdx(1) + 1 : cIdx(2), iTrials); %handle period
%     newVc(:, segFrames(2) + 1 : segFrames(2) + (diff(cIdx(2:3))), iTrials) = Vc(:, cIdx(2) + 1 : cIdx(3), iTrials); %stimulus period
%     newVc(:, segFrames(3) + 1 : segFrames(3) + (diff(cIdx(3:4))), iTrials) = Vc(:, cIdx(3) + 1 : cIdx(4), iTrials); %delay period
%     
%     if segFrames(4) + (diff(cIdx(4:5))) > segFrames(5)
%         newVc(:, segFrames(4) + 1 : segFrames(5), iTrials) = Vc(:, cIdx(4) + 1 : cIdx(4) + (segFrames(5) - segFrames(4)), iTrials); %response period
%     else
%         newVc(:, segFrames(4) + 1 : segFrames(4) + (diff(cIdx(4:5))), iTrials) = Vc(:, cIdx(4) + 1 : cIdx(5), iTrials); %response period
%     end
% end