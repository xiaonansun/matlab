function bhv = rateDisc_loadDetectionBhv(Animals, cPath, newRun)

if ~exist('cPath', 'var') || isempty(cPath)
    cPath = '\\grid-hs\churchland_nlsas_data\\data\Behavior_Simon\';
% cPath = '\\CHURCHLANDNAS\homes\DOMAIN=CSHL\smusall\Behavior_Simon\';
end

if strcmpi(Animals, 'EMX')
%     Animals = {'mSM80' 'mSM81' 'mSM82' 'mSM83' 'mSM84' 'mSM85' 'mSM86'};
    Animals = {'mSM80' 'mSM81' 'mSM82' 'mSM85' 'mSM86'};
elseif strcmpi(Animals, 'FezF')
    Animals = {'Fez7' 'Fez11' 'Fez13' 'Fez17' 'Fez18' 'Fez19'};
elseif strcmpi(Animals, 'Plexin')
    Animals = {'Plex05' 'Plex06' 'Plex07' 'Plex08'};
elseif strcmpi(Animals, 'CSP')
    Animals = {'CSP7' 'CSP8' 'CSP20' 'CSP24' 'CSP25'};
elseif strcmpi(Animals, 'Control')
    Animals = {'CTP3' 'CTP7'};
end
    
if ~exist('newRun', 'var')
    newRun = false;
end
if ~newRun
    try
        load([cPath 'rateDisc' filesep 'optoDetect_' Animals{:}], 'bhv')
    catch ME
        disp(ME.message);
        newRun = true;
        fprintf('Couldnt load processed bhv data. Loading raw files instead.\n')
    end
end
        
if newRun
    bhv = [];
    for iAnimals = 1:length(Animals)
        [~,cBhv] = rateDisc_optoStim(Animals{iAnimals},cPath , inf, 0.75);
        
        if ~isempty(cBhv)
            cBhv.AnimalID = ones(1, length(cBhv.Rewarded)) * iAnimals;
            if ~isempty(bhv)
                cBhv.SessionNr = cBhv.SessionNr + max(bhv.SessionNr);
            end
            cBhv = rateDisc_checkOptoPower(cBhv, Animals{iAnimals}); %check for missing power values
            if (sum(isnan(cBhv.optoPower(cBhv.optoDur > 0)))) > 0
                error
            end
            bhv = appendBehavior(bhv,cBhv); %append into larger array
        end
    end
    bhv = selectBehaviorTrials(bhv,~ismember(bhv.SessionNr, unique(bhv.SessionNr(bhv.DistStim > 0)))); %only use sessions that don't include discrmination trials
    if ~exist([cPath 'rateDisc' filesep], 'dir')
        mkdir([cPath 'rateDisc' filesep]);
    end
    
    %remove very large fields and save to file
    bhv = rmfield(bhv,'stimEvents');
    bhv = rmfield(bhv,'TrialSettings');
    save([cPath 'rateDisc' filesep 'optoDetect_' Animals{:}], 'bhv', '-v7.3')
end

%add animal names as output
bhv.Animals = Animals;
