% code to assess the beta kernels for all animals in a specific task
% period. This is meant to look at results during discirmination.
% Episode convention in recIdx: 
% 1 = early audio learning
% 2 = 5 to 50 percentile audio detection
% 3 = 50 to 95 percentile audio detection
% 4 = remaining audio detection
% 5 = Audio discrimination
% 6 = Novice tactile discrimination
% 7 = 5 to 50 percentile tactile detection
% 8 = 50 to 95 percentile tactile detection
% 9 = remaining tactile detection
% 10 = tactile discrimination
% 11 = everything else (usually mixed discriminnation)

[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
% cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
animals = dataOverview(:,1);
fileExt = 'org'; %non-orthogonalized model
load('allenDorsalMapSM.mat');
allenMask = dorsalMaps.allenMask;
load('trimMap.mat');
trimMap = arrayShrink(trimMap, allenMask, 'merge');
allenMask = dorsalMaps.allenMask;
areaIdx = [23, 33, 11, 5]; %index for left auditory, visAM, SShindlimb and M2 from trimMap
taskPeriod = 4; 
% cRegs = {'Choice' 'lAudStim' 'rAudStim' 'prevChoice' 'prevReward' 'nextChoice'}; %regressors of interest
cRegs = {'lfirstAudStim' 'rfirstAudStim'}; %regressors of interest
groups = {'mSM', 'Fez', 'Plex'};

%% get some data
betaData = cell(length(cRegs),length(groups));
groupCnt = zeros(1,length(groups));
for iAnimals = 1 : 4

    cAnimal = animals{iAnimals}; % current animal
    bPath = [cPath cAnimal filesep 'blockData' filesep]; % path for blockdata
    cGroup = NaN;
    for iGroups = 1 : length(groups)
        if contains(cAnimal, groups{iGroups})
            cGroup = iGroups; %animal is member of this group
        end
    end
    
    if ~isnan(cGroup)
        groupCnt(cGroup) = groupCnt(cGroup) + 1;

        % get beta kernels
        load([bPath fileExt 'betaMaps_episode1_reg1.mat'],'regLabels'); %get labels
        
        for iRegs = 1 : length(cRegs)
            cIdx = find(ismember(regLabels, cRegs{iRegs}));
            load([bPath fileExt 'betaMaps_episode' num2str(taskPeriod) '_reg' num2str(cIdx) '.mat'],'cBeta');
            
            if ~isempty(betaData{iRegs,cGroup}) && size(cBeta,2) ~= size(betaData{iRegs,cGroup},2)
                if size(cBeta,2) < size(betaData{iRegs,cGroup},2)
                    betaData{iRegs,cGroup} = betaData{iRegs,cGroup}(:, 1:size(cBeta,2));
                elseif size(cBeta,2) > size(betaData{iRegs,cGroup},2)
                    cBeta = cBeta(:, 1:size(betaData{iRegs,cGroup},2));
                end
            end
            betaData{iRegs,cGroup} = runMean(betaData{iRegs,cGroup}, cBeta, groupCnt(cGroup));
        end
        fprintf('Done - %s. %d / %d\n', cAnimal, iAnimals, length(animals));
    end
end

