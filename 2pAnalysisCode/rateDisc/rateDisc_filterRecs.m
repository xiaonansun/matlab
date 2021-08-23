function [recs, useIdx] = rateDisc_filterRecs(trainDates, recs, trainDur)
% code to selected recordings from a certain time range.

if ~isempty(trainDates)
    [cTrainDates, cTrainLabels] = rateDisc_convertTrainDates(trainDates);
else
    cTrainLabels = {};
end

if sum(ismember(cTrainLabels, trainDur)) >= 2
    if sum(ismember(cTrainLabels, trainDur)) > 2
        warning('!! More than 2 entries found for training range. Data will span entire range between earliest and latest date !!')
    end
    
    cDates = datenum(cTrainDates(contains(cTrainLabels, trainDur)));
    fprintf('Selected range: %s; First date: %s, Last date: %s.\n', trainDur, datestr(min(cDates)), datestr(max(cDates))); 

elseif strcmpi(trainDur, 'audiolearn')
    cDates = [-inf datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-50thPrc')))];
    fprintf('Selected range: %s; Using all dates until: %s.\n', trainDur, datestr(max(cDates)));
    
elseif strcmpi(trainDur, 'audioDetect')
    cDates = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-50thPrc'))) datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Last')))];
    fprintf('Selected range: %s; Using all dates from %s until: %s.\n', trainDur, datestr(min(cDates)), datestr(max(cDates)));
    
elseif strcmpi(trainDur, 'audioPerform')
    cDates = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-50thPrc'))) datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-last')))];
    fprintf('Selected range: %s; Using all dates from %s until: %s.\n', trainDur, datestr(min(cDates)), datestr(max(cDates)));
    
elseif strcmpi(trainDur, 'audioDisc')
    cDates = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-first')))  datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-last')))];
    fprintf('Selected range: %s; Using all dates from %s until: %s.\n', trainDur, datestr(min(cDates)), datestr(max(cDates)));
        
elseif strcmpi(trainDur, 'allAudio')
    cDates = [-inf datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-last')))];
    fprintf('Selected range: %s; Using all dates until: %s.\n', trainDur, datestr(max(cDates)));
    
elseif strcmpi(trainDur, 'allTac')
    cDates = [datenum(cTrainDates(contains(cTrainLabels, 'TacDetect-Learn-5thPrc'))) datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-last')))];
    fprintf('Selected range: %s; Using all dates until: %s.\n', trainDur, datestr(max(cDates))); 
    
elseif strcmpi(trainDur, 'rateDisc')
    cDates = [-inf max(datenum(cTrainDates))];
    fprintf('Selected range: %s; Using all recordings until last training date.\n', trainDur);
    
elseif strcmpi(trainDur, 'all')
    cDates = [-inf inf];
    fprintf('Selected range: %s; Using all available data.\n', trainDur);
else
    error('Unknown training range. Use "all" to get all data.')
end

useIdx = false(1, length(recs));
dateIdx = NaN(1, length(recs));
for iRecs = 1 : length(recs)
    try
        useIdx(iRecs) = datenum(recs(iRecs).name) >= min(cDates) && datenum(recs(iRecs).name) <= max(cDates);
        dateIdx(iRecs) = datenum(recs(iRecs).name);
    catch
        useIdx(iRecs) = false;
    end
end
dateIdx = dateIdx(useIdx);
recs = recs(useIdx);
[~, sortInd] = sort(dateIdx, 'ascend');
recs = recs(sortInd);
