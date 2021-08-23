%function rateDisc_alignUtoNAS
cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
tPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\BpodImager\Animals\'; %churchlandNAS path

load('allenDorsalMapSM.mat');
allenMask = dorsalMaps.allenMask;
[xRange, yRange] = rateDisc_maskRange(dorsalMaps.allenMask); % get inner range of allenMask

% get animal info
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates, allRegs] = rateDiscRecordings;
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
animals = dataOverview(:,1);
animals = animals(1:10);

for iAnimals = 1 : length(animals)

    % get recs
    cAnimal = animals{iAnimals};
    recs = dir([cPath cAnimal filesep 'SpatialDisc' filesep]);
    recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs, trainingRange); %this sorts recordings by date
%     recs = cat(1,recs(:).name);
%     save([tPath cAnimal filesep 'recs_' trainingRange '.mat'], 'recs')
    
    for iRecs = 1 : length(recs)
        % align U and save to churchlandNAS
        try
            vcFile = 'rsVc.mat';
            load([cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep vcFile], 'U', 'Vc')
        catch
            vcFile = 'Vc.mat';
            load([cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep vcFile], 'U', 'Vc')
        end
        load([cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep 'opts2.mat'], 'opts')
        
        U = alignAllenTransIm(single(U),opts.transParams); %align to allen
        U = single(rateDisc_removeOutline(U,10)); %remove potential artifact from alignment
        mask = ~isnan(U(:,:,1));
        
        Vc = reshape(Vc, size(Vc,1), []);
        nanIdx = isnan(mean(Vc, 1));
        [q, r] = qr(Vc(:,~nanIdx)', 0); %get qr results
        
        if ~exist([tPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep],'dir')
            mkdir([tPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]);
        end
        
        save([tPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep 'mask.mat'], 'mask');
        save([tPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep 'alignU.mat'], 'U', '-v7.3');
        save([tPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep 'QR.mat'], 'q', 'r', 'nanIdx', '-v7.3');
        copyfile([cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep vcFile], ...
            [tPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep vcFile]);
        
        fprintf('Current recording: %d-%d \n', iRecs, length(recs));
    end
end