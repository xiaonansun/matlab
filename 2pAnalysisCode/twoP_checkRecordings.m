
[dataOverview, ~, ~, ~, ~, ~, ~, cPath] = twoP_delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);


for iAnimals = 70 : length(animals)
    
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep]; disp(fPath); %current data path
    load([fPath 'data.mat']);
    
    rejIdx = twoP_checkTrialDrift(data, dataOverview{iAnimals, 9}(1), dataOverview{iAnimals, 9}(2), dataOverview{iAnimals, 6}); %check for drift over trials
    data.rejIdx = rejIdx;
    
    if ~any(isinf(dataOverview{iAnimals, 6}))
        bTrials(~ismember(bTrials,data.trialNumbers(dataOverview{iAnimals, 6}))) = []; %reject trials based on pre-determiend trial index. This is in case something bad happened in the session.
    end
    data.bhvTrials = bTrials;

    save([fPath 'data.mat'], 'data', 'bTrials');
    close all; clear data;

    twoP_prepareData(fPath,dataOverview{iAnimals,4},dataOverview{iAnimals,6}); %prepare 2p data and trials and move to target folder
%     delayDec_RegressModel(cPath,animals{iAnimals},recs{iAnimals},'twoP');
    
%     delayDec_testRegs(cPath, animals{iAnimals}, recs{iAnimals}, true, 'twoP')
%     delayDec_testRegs(cPath, animals{iAnimals}, recs{iAnimals}, false, 'twoP')
% 
%     disp('BodyVars'); Behavior_computeBodyVars(cPath, animals{iAnimals}, recs{iAnimals}, 2)
%     disp('PupileVars'); Behavior_computePupilVars([fPath 'BehaviorVideo'],false)
%     twoP_delayRegressModel(cPath,animals{iAnimals},recs{iAnimals});
%     twoP_testRegs(cPath,animals{iAnimals},recs{iAnimals},true)
%     twoP_testRegs(cPath,animals{iAnimals},recs{iAnimals},false)
    
%     Behavior_redoCombinedSVD([fPath 'BehaviorVideo\'])
%     Behavior_computeFaceVars(cPath, animals{iAnimals}, recs{iAnimals}, 1, false);
%     Behavior_computePupilVars([fPath 'BehaviorVideo'],true,true)
end