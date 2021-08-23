function selTrials = Behavior_AnalysisTrialselection(SessionData,trials,mod1,mod2,rewardType)
% code to select trials that are balanced, as required for later analysis
% Currently balances trialcounts for left vs right and mod1 vs mod2

if ~exist('rewardType','var')
    rewardType = true; %only use correct trials by default
end

%% get reaction times
rTimes = Behavior_ReactionTimes(SessionData); %compute reaction times
rTimes = rTimes(trials);

%% get modality indices
sucInd = SessionData.Rewarded(trials) == rewardType & SessionData.Assisted(trials) ...  %find succesful unisensory trials
    & (SessionData.TouchCnt(trials) < 6); % reject trials with more than 6 baseline lever touches
modInd{1} = SessionData.StimType(trials) == mod1;
modInd{2} = SessionData.StimType(trials) == mod2;

%% find even amount of trials for both sides / modalities
selTrials = []; %trials to be kept
for iSides = 1:2
    
   ind{1} = SessionData.ResponseSide(trials) == iSides & modInd{1} & sucInd; %correct trials for a given side/modality
   ind{2} = SessionData.ResponseSide(trials) == iSides & modInd{2} & sucInd;
   [tCnt(iSides), minMod(iSides)] = min([sum(ind{1}) sum(ind{2})]);
   
end
[~,b] = min(tCnt); %minimum amount of trials for each side/modality

% get side/modality which has the least responses and match all other cases to its response time distribution
compResp = rTimes(SessionData.ResponseSide(trials) == b & SessionData.StimType(trials) == minMod(b) & sucInd);

for iSides = 1:2
    for iMods = 1:2
        for iResp = 1 : length(compResp) %find response times from higher trialcount population that match lower trialcounts best

            ind = find(SessionData.ResponseSide(trials) == iSides & modInd{iMods} & sucInd); %correct trials for a given side/modality
            tempInd = true(1,length(ind)); %all trials for current case
            tempInd(ismember(ind,selTrials)) = false; %dont re-use trials that were used already
            ind = ind(tempInd); %index of trials that havent been selected yet
            
            [~,tempResp] = min(abs(rTimes(ind) - compResp(iResp))); %find closest response time
            selTrials = [selTrials ind(tempResp)];
            
        end
    end
end
