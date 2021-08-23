[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, ~, fPath, trainDates] = rateDiscRecordings;
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
animals = dataOverview(:,1);
recs = dataOverview(:,3);

iAnimals = 1;
animal = animals{iAnimals};

recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, 'rateDisc'); %this sorts recordings by date
