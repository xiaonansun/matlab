function rateDisc_batchRunModel(animal)

if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
end

[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, 'rateDisc'); %this sorts recordings by date

for iRecs = 1 : size(recs,1)
    
%     fPath = [cPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
    
%     if exist([fPath 'interpVc.mat'],'file') == 2
%         load([fPath 'interpVc.mat'],'frames');
%     else
%         frames = 0;
%     end
%     
%     if frames ~= 150
        disp([animal ' - ' recs(iRecs).name]);
        
        cLine = ['qsub -l m_mem_free=2G -pe threads 8 -binding linear:8 runModel.sh ' ...
            animal ' ' num2str(iRecs)];
        system(cLine);
%     end
end