cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
% cPath = '\\CHURCHLANDNAS\homes\DOMAIN=CSHL\smusall\Behavior_Simon\';
% Animals = {'mSM83' 'mSM84'};
% Animals = {'mSM80' 'mSM81' 'mSM82'};
% Animals = {'Plex05' 'Plex06'};
Animals = {'Fez18' 'Fez19'};
% Animals = {'CSP7' 'CSP8' 'CSP20' 'CSP24' 'CSP25'};
% Animals = {'Fez7' 'Fez11' 'Fez13' 'Fez17' 'Fez18' 'Fez19'};
% 
modLabels = {'all'}; %labels for different modalities
timeLabels = {'EarlyStim' 'Delay'}; %labels for different opto stim times
siteLabels = {'parietal' 'frontal'}; %labels for different opto locations
modId = 2; %id for modalities. Default is audio, somatosensory and audiosomatosensory
stimSides = 3;
distBins = 6;
newRun = false;

if ~newRun
    try
        load([cPath 'rateDisc' filesep 'optoBhv_Discrimination_' Animals{:}], 'bhv')
    catch ME
        disp(ME.message);
        newRun = true;
        fprintf('Couldnt load processed bhv data. Loading raw files instead.\n')
    end
end
        
if newRun
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
    
    if ~exist([cPath 'rateDisc' filesep], 'dir')
        mkdir([cPath 'rateDisc' filesep]);
    end
    
    %remove very large fields and save to file
    bhv = rmfield(bhv,'TrialSettings');
    save([cPath 'rateDisc' filesep 'optoBhv_Discrimination_' Animals{:}], 'bhv', '-v7.3')
end

%% add light power information
lowPowerOn = '18-Jun-2018'; %date when light power was reduced
% lowPowerOn = '31-Dec-2019'; %date when light power was reduced
bhv.optoPower = repmat(10,1,length(bhv.Rewarded)); %default power is 10mW for 0.4mm implant
bhv.optoPower(bhv.date > datenum(lowPowerOn)) = 5; %reduced power is only 5 mW

stimLoc = [1 2]; %index for locations, use two values here. 1 = frontal, 2 = parietal, 3 = S1

%% compute discrimination performance - Bilateral stim and delay stimuation
oInd = bhv.optoDur == 1.5 & bhv.optoSide == 3 & bhv.optoPower == max(bhv.optoPower); %bilateral stim and delay stimuation at max. power.
pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & ismember(bhv.AnimalID,2); %only use active trials

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





%% compute left/right discrimination performance - Unilateral stim and delay stimuation at high light power
figure
for iAnimals = 1 : length(Animals) + 1
    
    if iAnimals > length(Animals)
        oInd = bhv.optoDur <= 1.5 &  bhv.optoDur > 0 & bhv.optoSide ~= 3 & bhv.optoPower == max(bhv.optoPower); %bilateral stim and delay stimuation at max. power.
    else
        oInd = bhv.optoDur <= 1.5 &  bhv.optoDur > 0 & bhv.optoSide ~= 3 & bhv.optoPower == max(bhv.optoPower) & bhv.AnimalID == iAnimals; %bilateral stim and delay stimuation at max. power.
    end
    
    pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %only use active trials
    
    clear cInd
    cInd(2,:) = pInd & oInd & bhv.stimLocation == stimLoc(1) & bhv.optoType == 1 & bhv.optoSide == 1; %optogenetic; left frontal stim only
    cInd(3,:) = pInd & oInd & bhv.stimLocation == stimLoc(1) & bhv.optoType == 1 & bhv.optoSide == 2; %optogenetic; right frontal stim only
    cInd(1,:) = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(cInd(2,:)))); %non-optogenetic trials from same sessions
    
    cInd(5,:) = pInd & oInd & bhv.stimLocation == stimLoc(2) & bhv.optoType == 1 & bhv.optoSide == 1; %optogenetic; left parietal stim
    cInd(6,:) = pInd & oInd & bhv.stimLocation == stimLoc(2) & bhv.optoType == 1 & bhv.optoSide == 2; %optogenetic; right parietal stim
    cInd(4,:) = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(cInd(5,:)))); %non-optogenetic trials from same sessions
    
    subplot(length(Animals)+1,2,(iAnimals-1)*2 + 1);
    cColors = {'k' 'g' 'r' 'k' 'g' 'r'}; clear h
    Cnt = 0;
    for x = 1 : 3
        
        Cnt = Cnt +1;
        cInd(x,:) = rateDisc_equalizeTrialsPerMouse(bhv, cInd(x,:)); %equalize trial counts for animals in current selection
        [distRatio, rightChoice, nTrials, params, cFit] = rateDisc_audioDiscCurve(bhv, cInd(x,:), distBins);
        
        %Get 95% Wilson binomial CIs for data and plot performance
        z = 1.96;
        pChoseHigh = rightChoice./nTrials;
        dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        
        %Plot fit and real data
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
        h(Cnt) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
        
    end
    legend(h,{'Control','Left','Right frontal Stim+Delay'},'location','northwest');
    xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
    xlabel('Distractor ratio');
    ylim([0 1]); ylabel('Proportion chose right')
    hold off
    axis square
    if iAnimals > length(Animals)
         title('Frontal Stim+Delay - All mice');
    else
         title(['Frontal Stim+Delay - ' Animals{iAnimals}]);
    end
    
    subplot(length(Animals)+1,2,(iAnimals-1)*2 + 2); Cnt = 0;
    for x = 4 : 6
        
        Cnt = Cnt +1;
        cInd(x,:) = rateDisc_equalizeTrialsPerMouse(bhv, cInd(x,:)); %equalize trial counts for animals in current selection
        [distRatio, rightChoice, nTrials, params, cFit] = rateDisc_audioDiscCurve(bhv, cInd(x,:), distBins);
        
        %Get 95% Wilson binomial CIs for data and plot performance
        z = 1.96;
        pChoseHigh = rightChoice./nTrials;
        dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        
        %Plot fit and real data
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
        h(Cnt) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
        
    end
%     legend(h,{'Control','Left parietal Stim+Delay','Right parietal Stim+Delay'},'location','northwest');
    xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
    xlabel('Distractor ratio');
    ylim([0 1]); ylabel('Proportion chose right')
    hold off
    axis square
    if iAnimals > length(Animals)
         title('Parietal Stim+Delay - All mice');
    else
         title(['Parietal Stim+Delay - ' Animals{iAnimals}]);
    end
end

%% compute stim/delay discrimination performance - Bilateral stimulation at high power
figure
for iAnimals = 1 : length(Animals) + 1
    
    if iAnimals > length(Animals)
        oInd = bhv.optoDur < 1.5 & bhv.optoSide == 3 & bhv.optoPower == min(bhv.optoPower); %bilateral stim and delay stimuation at max. power.
    else
        oInd = bhv.optoDur < 1.5 & bhv.optoSide == 3 & bhv.optoPower == min(bhv.optoPower) & bhv.AnimalID == iAnimals; %bilateral stim and delay stimuation at max. power.
    end
    pInd = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %only use active trials
    
    clear cInd
    cInd(2,:) = pInd & oInd & bhv.stimLocation == stimLoc(1) & bhv.optoType == 1; %optogenetic; stimulus - frontal
    cInd(3,:) = pInd & oInd & bhv.stimLocation == stimLoc(1) & bhv.optoType == 2; %optogenetic; delay - frontal
    cInd(1,:) = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(cInd(2,:)))); %non-optogenetic trials from same sessions
    
    cInd(5,:) = pInd & oInd & bhv.stimLocation == stimLoc(2) & bhv.optoType == 1; %optogenetic; stimulus - parietal stim
    cInd(6,:) = pInd & oInd & bhv.stimLocation == stimLoc(2) & bhv.optoType == 2; %optogenetic; delay - parietal stim
    cInd(4,:) = pInd & bhv.optoDur == 0 & ismember(bhv.SessionNr, unique(bhv.SessionNr(cInd(5,:)))); %non-optogenetic trials from same sessions
    
    subplot(length(Animals)+1,2,(iAnimals-1)*2 + 1);
%     subplot(1,2,1);
    cColors = {'k' 'g' 'r' 'k' 'g' 'r'}; clear h
    Cnt = 0;
    for x = 1 : 3
        
        Cnt = Cnt +1;
        cInd(x,:) = rateDisc_equalizeTrialsPerMouse(bhv, cInd(x,:)); %equalize trial counts for animals in current selection
        [distRatio, rightChoice, nTrials, params, cFit] = rateDisc_audioDiscCurve(bhv, cInd(x,:), distBins);
        
        %Get 95% Wilson binomial CIs for data and plot performance
        z = 1.96;
        pChoseHigh = rightChoice./nTrials;
        dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        
        %Plot fit and real data
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
        h(Cnt) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
        
    end
%     legend(h,{'Control','Stim frontal','Delay frontal'},'location','northwest');
    xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
    xlabel('Distractor ratio');
    ylim([0 1]); ylabel('Proportion chose right')
    hold off
    axis square
    if iAnimals > length(Animals)
         title('Frontal Stim / Delay - All mice');
    else
         title(['Frontal Stim / Delay - ' Animals{iAnimals}]);
    end
    
    subplot(length(Animals)+1,2,(iAnimals-1)*2 + 2); Cnt = 0;
%     subplot(1,2,2);
    for x = 4 : 6
        
        Cnt = Cnt +1;
        cInd(x,:) = rateDisc_equalizeTrialsPerMouse(bhv, cInd(x,:)); %equalize trial counts for animals in current selection
        [distRatio, rightChoice, nTrials, params, cFit] = rateDisc_audioDiscCurve(bhv, cInd(x,:), distBins);
        
        %Get 95% Wilson binomial CIs for data and plot performance
        z = 1.96;
        pChoseHigh = rightChoice./nTrials;
        dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
        
        %Plot fit and real data
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', cColors{x}); hold on
        h(Cnt) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', cColors{x}, 'MarkerFaceColor','w','linewidth',2);
        
    end
    legend(h,{'Control','Stimulus','Delay'},'location','northwest');
    xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
    xlabel('Distractor ratio');
    ylim([0 1]); ylabel('Proportion chose right')
    hold off
    axis square
    if iAnimals > length(Animals)
         title('Parietal Stim / Delay - All mice');
    else
         title(['Parietal Stim / Delay - ' Animals{iAnimals}]);
    end
end