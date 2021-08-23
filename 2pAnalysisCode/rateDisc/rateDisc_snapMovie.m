% rateDisc_snapMovie

[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, ~, fPath, trainDates] = rateDiscRecordings;
cPath = 'Y:\data\BpodImager\Animals\';
stepSize = 15;
regType = 'ridge'; %lasso or ridge
% animals = dataOverview(:,1);
animals = {'mSM63';'mSM64';'mSM65';'mSM66';'Fez7';'Fez10';'Plex01';'Plex02'};
highThresh = 1; %threshold for colormap to increase contrast
winDiam = 200; %diameter of circle for second alignment

for iAnimals = 5 : length(animals)
    %% get recordings and sort by date
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep];
    recs = ls(fPath);
    recs = recs(~ismember(recs(:,1), '.'), :);
    
    load([fPath strtrim(recs(2,:)) filesep 'Snapshot_1.mat'])
    snapMovie = NaN(size(snap,1), size(snap,2), size(recs,1), 'single');
    
    cIdx = false(1, size(recs,1));
    for iRecs = 1 : size(recs,1)
        try            
            load([fPath strtrim(recs(iRecs,:)) filesep 'Snapshot_1.mat']);
            snapMovie(:, :, iRecs) = snap;
        catch cIdx(iRecs) = true;
        end
    end
    recs(cIdx, :) = [];
    snapMovie(:, :, cIdx) = [];
    clear cIdx
    
    %% save as mat and h5 files
    if exist([cPath animals{iAnimals} filesep 'snapMovie.h5\'], 'file')
        delete([cPath animals{iAnimals} filesep 'snapMovie.h5\']);
    end
    
    recIdx = ~isnan(snapMovie(1,1,:)); %check for omitted recordings
    save([cPath animals{iAnimals} filesep 'snapMovie.mat'], 'snapMovie', 'recIdx');
    
    snapMovie = arrayResize(snapMovie, 2);
    snapMovie(272, :, :) = 0; %this will keep leap from crashing :)
    
    snapMovie = snapMovie ./ (2^16 - 1);
    snapMovie(snapMovie > highThresh) = highThresh;
    snapMovie = reshape(snapMovie(:, :, recIdx), size(snapMovie,1), size(snapMovie,2), 1, sum(recIdx));
    snapMovie = uint8(snapMovie .* 255);
    
    h5create([cPath animals{iAnimals} filesep 'snapMovie.h5\'],'/box', size(snapMovie),'ChunkSize', size(snapMovie),'Datatype','uint8')
    h5write([cPath animals{iAnimals} filesep 'snapMovie.h5\'], '/box', snapMovie)
    
    fprintf('Done: %s\n', animals{iAnimals});
    
end

%% combine all snaps into larger stack and save separately
allSnap = cell(1, length(animals));
for iAnimals = 1 : length(animals)
    load([cPath animals{iAnimals} filesep 'snapMovie.mat'], 'snapMovie', 'recIdx');
    allSnap{iAnimals} = snapMovie;
end

allSnap = cat(3, allSnap{:});
save([cPath filesep 'allSnap.mat'], 'allSnap');

allSnap = arrayResize(allSnap, 2);
allSnap(272, :, :) = 0; %this will keep leap from crashing :)

allSnap = allSnap ./ (2^16 - 1);
allSnap(allSnap > highThresh) = highThresh;
allSnap = reshape(allSnap, size(allSnap,1), size(allSnap,2), 1, size(allSnap,3));
allSnap = uint8(allSnap .* 255);

h5create([cPath filesep 'allSnap.h5\'],'/box', size(allSnap),'ChunkSize', size(allSnap),'Datatype','uint8')
h5write([cPath filesep 'allSnap.h5\'], '/box', allSnap)

