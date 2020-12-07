% function twoP_batchRun

[dataOverview, ~, ~, ~, ~, ~, ~, cPath] = twoP_delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);


for iAnimals = 1 : length(animals)
% for iAnimals = [5 6 34 35 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62]
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep]; disp(fPath); %current data path
%     twoP_prepareData(fPath,dataOverview{iAnimals,4},dataOverview{iAnimals,6}); %prepare 2p data and trials and move to target folder
    
try
    load([fPath 'interpVc.mat'],'DS','frames','Vc') %Vc that was used for the model
    [b, a] = butter(2, 0.01/31, 'high');
    DS = filtfilt(b, a, DS');
    DS = bsxfun(@minus, DS, median(DS));    
    moveIdx = any(DS' > 2); %don't use data with motion above 1 standard deviation
    
    dataCheck(iAnimals, 1) = length(moveIdx);
    dataCheck(iAnimals, 2) = sum(moveIdx);
    disLoc(iAnimals,:) = nanstd(DS, [], 1);
    
    dataSize(iAnimals, 1) = size(Vc,2) / frames;
    dataSize(iAnimals, 2) = size(Vc,1);
    
%     delayDec_RegressModel(cPath,animals{iAnimals},recs{iAnimals},'twoP');
%     close all;
catch
    disp('Current recording failed');
end
end

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
% end