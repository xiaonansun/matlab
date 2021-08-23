cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
% Animals = {'mSM83' 'mSM84'};
% Animals = {'mSM80' 'mSM81' 'mSM82'};
% Animals = {'Plex05' 'Plex06'};
Animals = {'mSM80'};
% Animals = {'Fez11'};
modLabels = {'all'}; %labels for different modalities
timeLabels = {'EarlyStim' 'Delay'}; %labels for different opto stim times
siteLabels = {'parietal' 'frontal'}; %labels for different opto locations
modId = 2; %id for modalities. Default is audio, somatosensory and audiosomatosensory
stimSides = 3;
distBins = 6;

bhv = [];
for iAnimals = 1 : length(Animals)
    [Performance,cBhv] = rateDisc_optoStim(Animals{iAnimals},cPath , inf, 0.6, true);
    
    if ~isempty(cBhv)
        cBhv.AnimalID = ones(1, length(cBhv.Rewarded)) * iAnimals;
        if ~isempty(bhv)
            cBhv.SessionNr = cBhv.SessionNr + max(bhv.SessionNr);
        end
        bhv = appendBehavior(bhv,cBhv); %append into larger array
    end
end
bhv = selectBehaviorTrials(bhv,ismember(bhv.SessionNr, unique(bhv.SessionNr(bhv.DistStim > 0)))); %only use sessions that include discrmination trials

% add light power information
lowPowerOn = '18-Jun-2019'; %date when light power was reduced
% lowPowerOn = '31-Dec-2019'; %date when light power was reduced
bhv.optoPower = repmat(10,1,length(bhv.Rewarded)); %default power is 10mW for 0.4mm implant
bhv.optoPower(bhv.date > datenum(lowPowerOn)) = 5; %reduced power is only 5 mW

stimLoc = [1 2]; %index for locations, use two values here. 1 = frontal, 2 = parietal, 3 = S1

%% compute discrimination performance - Bilateral stim and delay stimuation

oInd = bhv.optoDur == 1.5 & bhv.optoSide == 3 & bhv.optoPower == max(bhv.optoPower); %bilateral stim and delay stimuation at max. power.
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %only use active trials

clear cInd
cInd(2,:) = pInd & oInd & bhv.stimLocation == stimLoc(1) & bhv.optoType == 1; %optogenetic; frontal stim only
cInd(1,:) = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(cInd(2,:)))); %non-optogenetic trials from same sessions
cInd(4,:) = pInd & oInd & bhv.stimLocation == stimLoc(2) & bhv.optoType == 1; %optogenetic; parietal stim only
cInd(3,:) = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(cInd(4,:)))); %non-optogenetic trials from same sessions
% cInd(5,:) = pInd & oInd & bhv.stimLocation == stimLoc(2) & bhv.optoType == 2; %optogenetic; parietal delay only

figure
subplot(1,2,1);
cColors = {'k' 'b' 'k' 'r'}; clear h
for x = 1 : 2
    
    cInd(x,:) = rateDisc_equalizeTrialsPerMouse(bhv, cInd(x,:)); %equalize trial counts for animals in current selection
    [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, cInd(x,:), distBins);
    
    %Plot fit and real data
    plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
    h(x) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
    
end
legend(h,{'Control','Frontal Stim+Delay'},'location','northwest');
xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
ylim([0 1]); ylabel('Proportion chose right'); hold off; axis square
xlabel('Distractor ratio'); title('Frontal bilateral inactivation');
disp(nTrials);

subplot(1,2,2);
for x = [3 4]
    
    cInd(x,:) = rateDisc_equalizeTrialsPerMouse(bhv, cInd(x,:)); %equalize trial counts for animals in current selection
    [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(bhv, cInd(x,:), distBins);
    
    %Plot fit and real data
    plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
    h(x) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
    
end
legend(h([3 4]),{'Control','Parietal Stim+Delay'},'location','northwest');
xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
ylim([0 1]); ylabel('Proportion chose right'); hold off; axis square
xlabel('Distractor ratio'); title('Parietal bilateral inactivation');
disp(nTrials);