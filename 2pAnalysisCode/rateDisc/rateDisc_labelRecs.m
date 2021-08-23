function [recIdx,dateIdx] = rateDisc_labelRecs(trainDates, recs)
% code to create an idex that identifies which phase of training each
% recording belongs to. 
% Convention in recIdx: 
% 1 = early audio learning
% 2 = 5 to 50 percentile audio detection
% 3 = 50 to 95 percentile audio detection
% 4 = remaining audio detection
% 5 = Audio discrimination
% 6 = Novice tactile discrimination
% 7 = 5 to 50 percentile tactile detection
% 8 = 50 to 95 percentile tactile detection
% 9 = remaining tactile detection
% 10 = tactile discrimination
% 11 = everything else (usually mixed discriminnation)

recIdx = zeros(1,length(recs)); %index for each recording.

if ~isempty(trainDates)
    [cTrainDates, cTrainLabels] = rateDisc_convertTrainDates(trainDates);
else
    cTrainLabels = {};
end

% get index for different training ranges
dIdx(1,:) = [-inf datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-5thPrc')))-1];                                                                          % 1 = 5 to 50 percentile audio detection
dIdx(2,:) = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-5thPrc'))) datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-50thPrc')))];       % 2 = 5 to 50 percentile audio detection
dIdx(3,:) = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-50thPrc')))+1 datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-95thPrc')))];    % 3 = 50 to 95 percentile audio detection
dIdx(4,:) = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDetect-Learn-95thPrc')))+1 datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-first')))-1];            % 4 = remaining audio detection
dIdx(5,:) = [datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-first'))) datenum(cTrainDates(contains(cTrainLabels, 'AudioDisc-last')))];                           % 5 = remaining audio detection
try
    dIdx(6,:) = [datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-Novice-first'))) datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-Novice-last')))];                 % 6 = Novice tactile discrimination
    dIdx(7,:) = [datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-Novice-last')))+1 datenum(cTrainDates(contains(cTrainLabels, 'TacDetect-Learn-50thPrc')))];           % 7 = 5 to 50 percentile tactile detection
    dIdx(8,:) = [datenum(cTrainDates(contains(cTrainLabels, 'TacDetect-Learn-50thPrc')))+1 datenum(cTrainDates(contains(cTrainLabels, 'TacDetect-Learn-95thPrc')))];        % 8 = 50 to 95 percentile tactile detection
    dIdx(9,:) = [datenum(cTrainDates(contains(cTrainLabels, 'TacDetect-Learn-95thPrc')))+1 datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-first')))-1];                % 9 = 50 to 95 percentile tactile detection
    dIdx(10,:) = [datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-first'))) datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-last')))];                              % 10 = tactile discrimination
    dIdx(11,:) = [datenum(cTrainDates(contains(cTrainLabels, 'TacDisc-last')))+1 inf];                                                                                      % 11 = everything else (usually mixed discriminnation)
end

dateIdx = NaN(1, length(recs)); %this is too keep the datenumber for each recording. Useful for sorting ect.
for iRecs = 1 : length(recs)
    
    for iRange = 1 : size(dIdx,1)
        if (datenum(recs(iRecs).name) >= dIdx(iRange,1)) && (datenum(recs(iRecs).name) <= dIdx(iRange,2))
            recIdx(iRecs) = iRange;
            dateIdx(iRecs) = datenum(recs(iRecs).name);
        end
    end
end