[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, ~, fPath, trainDates] = rateDiscRecordings;
cPath = 'Y:\data\BpodImager\Animals\';
stepSize = 15;
regType = 'ridge'; %lasso or ridge
animals = dataOverview(:,1);
recs = dataOverview(:,3);

% Cnt = 0; brokenRec = [];
for iAnimals = 1 : length(animals)
    %% get recordings and sort by date
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep];
    recs = ls(fPath);
    recs = recs(~ismember(recs(:,1), '.'), :);
    cDate = datenum(recs(:, 1:11)); %get dates from recordings
   
    if size(recs,2) > 11 %of there are multiple recordings from the same date
        for iRecs = find(ismember(recs(:,12), '_'))
           cDate(iRecs) = cDate(iRecs) + (str2double(recs(iRecs,13))/100);
        end
    end
    [cDate,ind] = sort(cDate,'ascend'); % sort by oldest first
    recs = recs(ind,:); %adjust order of filenames to get it to be chronological
    recs = mat2cell(recs, ones(1,size(recs,1)), size(recs,2)); %convert to cells
        
    for iRecs = 1:length(recs)
        rateDisc_RegressModel(cPath, animals{iAnimals}, recs{iRecs}, 'Widefield');
    end
end
        