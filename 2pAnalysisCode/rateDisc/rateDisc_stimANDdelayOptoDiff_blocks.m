function [optDetect, ctrlDetect, optError, ctrlError] = rateDisc_stimANDdelayOptoDiff_blocks(bhv,blockCnt)
%short function to compute detection performance with optogenetic
%inactiavtion during stim+delay period and compare performance against
%non-inactivation trials. Returns performance for each session and all
%trials combined.
% bhv is the behavioral array from bpod, blockCnt is the number of trials, used in each block.

%% for each animal, determine the early and late stage in training and get behavior for early or late inactivation
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.DistStim == 0; %use active detection trials on both sides
% oInd = (bhv.optoDur == 1.5 | bhv.optoDur == 1.7) & bhv.optoType == 1 & bhv.optoSide == 3 & bhv.DistStim == 0 & bhv.optoPower > 1; %bilateral stim+delay stimulation.
oInd = (bhv.optoDur == 0.5 | bhv.optoDur == 0.3) & bhv.optoSide == 3 & bhv.DistStim == 0 & bhv.optoPower == 10; %bilateral stim+delay stimulation.

for iAnimals = unique(bhv.AnimalID)
    aInd = pInd & ismember(bhv.AnimalID,iAnimals); %performed trials for current animal

    % get performance for early and late inacivation at different locations. 1 = frontal, 2 = parietal
    for stimLoc = 1 : 2
        cInd = find(aInd & oInd & bhv.stimLocation == stimLoc); %optogenetic trials for current location and animal
        iBlocks = 1 : blockCnt : length(cInd); %index for blocks
        
        for x = 1 : length(iBlocks)
            
            if x == length(iBlocks)
                optTrials = cInd(end - blockCnt + 1 :end); %make sure last block has enough trials
            else
                optTrials = cInd(iBlocks(x) : iBlocks(x) + blockCnt - 1); %optogenetic trials in current block
            end
            ctrlTrials = find(ismember(bhv.SessionNr, unique(bhv.SessionNr(optTrials))) & ~oInd); %non-optogenetic trials from same sessions
            
            [optError{iAnimals,stimLoc}(x,1), optError{iAnimals,stimLoc}(x,2), optDetect{iAnimals,stimLoc}(x)] = ...
                    Behavior_wilsonError(sum(bhv.Rewarded(optTrials)), length(optTrials)); %performance with error
                
            [ctrlError{iAnimals,stimLoc}(x,1), ctrlError{iAnimals,stimLoc}(x,2), ctrlDetect{iAnimals,stimLoc}(x)] = ...
                    Behavior_wilsonError(sum(bhv.Rewarded(ctrlTrials)), length(ctrlTrials)); %performance with error            
        end
    end
end
          
    