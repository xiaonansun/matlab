function Behavior_checkCamData

sPath = 'U:\space_managed_data\BpodImager\Animals\'; %sourcepath
dataOverview = delayDecRecordings;
animals = dataOverview(:,1);

for iAnimals = 1 : size(dataOverview,1)

    fPath = [sPath dataOverview{iAnimals,1} filesep 'SpatialDisc' filesep dataOverview{iAnimals,3} filesep 'BehaviorVideo' filesep]; %current source datapath
    movieFiles1 = dir([fPath '*0001_1.mj2']);
    movieFiles2 = dir([fPath '*0001_2.mj2']);
    
    cFile = [fPath movieFiles1.name];
    v = VideoReader(cFile);
    bFrame(:,:,iAnimals,1) = single(arrayResize(readFrame(v),2)); clear v
    
    cFile = [fPath movieFiles2.name];
    v = VideoReader(cFile);
    bFrame(:,:,iAnimals,2) = single(arrayResize(readFrame(v),2)); clear v
    
    
end
