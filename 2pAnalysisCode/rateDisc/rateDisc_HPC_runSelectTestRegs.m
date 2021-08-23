function rateDisc_HPC_runSelectTestRegs
% code to run the 'Widefield_testRegs' command over all animals in the
% current dataset. Execute from bnbdev1. Needs updated 'delayDecRecordings' to get all
% recordings.

if ispc
    %     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
else
    %     cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/'; %path to grid
end
trainingRange = 'audioDisc'; %use this to run analysis only in a certain range of training

%% get animal info
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates, cRegs] = rateDiscRecordings;
animals = dataOverview(:,1);

for iAnimals = 9:10
    
    % get recs
    cAnimal = animals{iAnimals};
    recs = dir([cPath cAnimal filesep 'SpatialDisc' filesep]);
    recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs, trainingRange); %this sorts recordings by date
    
    
    for iRecs = 1 : length(recs)
        disp([cAnimal ' - ' recs(iRecs).name]);
        
        for iRegs = 1 : length(cRegs)
%             cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 testRegs.sh ' ...
%                 cAnimal ' ' recs(iRecs).name ' nomean ' cRegs{iRegs}];
%             system(cLine);
            
            cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 testRegs.sh ' ...
                cAnimal ' ' recs(iRecs).name ' org ' cRegs{iRegs}];
            system(cLine);
        end
    end
end