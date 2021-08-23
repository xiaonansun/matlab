function recs = rateDisc_getRecsForAnimal(animal, trainingRange, cPath)

if ~exist('trainingRange','var') || isempty(trainingRange)
    trainingRange = 'all';
end
if ~exist('cPath','var') || isempty(cPath)
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
end


[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, trainingRange); %this sorts recordings by date

end