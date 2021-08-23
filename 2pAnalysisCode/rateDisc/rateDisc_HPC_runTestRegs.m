function rateDisc_HPC_runTestRegs
% code to run the 'Widefield_testRegs' command over all animals in the
% current dataset. Execute from bnbdev1. Needs updated 'rateDiscRecordings' to get all
% recordings.

if ispc
%     cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
else
%     cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/'; %path to grid
end
trainingRange = 'allAudio';
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
animals = dataOverview(:,1);

%% get animal info
for iAnimals = 9 : 10
        
    % get recs
    cAnimal = animals{iAnimals};
%     try
%         load([cPath cAnimal filesep 'blockData' filesep 'trialInfo.mat'],'recs')
%     catch
%         disp('Couldnt find trialInfo in blockData folder. Moving all recs instead.');
        recs = dir([cPath cAnimal filesep 'SpatialDisc' filesep]);
        recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs, trainingRange); %this sorts recordings by date
%     end

    for iRecs = 1 : length(recs)
        disp([cAnimal ' - ' recs(iRecs).name]);
        
        cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 runModel.sh ' ...
            cAnimal ' ' num2str(iRecs) ' ' trainingRange];
        system(cLine);
        
%         cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 testRegs.sh ' ...
%             cAnimal ' ' recs(iRecs).name ' true'];
%         system(cLine);
        
%         cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 testRegs.sh ' ...
%             cAnimal ' ' recs(iRecs).name ' false'];
%         system(cLine);
        

    end
end