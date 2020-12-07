function ind = findInM(M, theDate, tetrode, cluster)
% ind = findInM(M, theDate, tetrode, cluster)
%
% Find an element of M corresponding to the specified unit. M is the M
% struct to search in. theDate may be any substring of the all_data
% filename (including a date such as '2012-04-02' or an entire filename
% including the '.mat' at the end). tetrode should be the tetrode number,
% and cluster should be the cluster number. ind will be the index in M of
% the matching entry. If no matches, ind will be empty. If multiple matches
% (because, for example, an incomplete theDate string was used), a warning
% will be issued and ind will have length >1.


if ~ischar(theDate) || ~isnumeric(tetrode) || ~isnumeric(cluster)
  error('findInM:badArgType', 'theDate must be a string, tetrode and cluster must be numbers');
end

MFiles = {M.filename};
MTetrodes = [M.tetrode];
MClusters = [M.cluster];

% Allow partial matches for date
dateMatchCell = strfind(MFiles, theDate);
dateMatch = ~cellfun(@isempty, dateMatchCell);

tetrodeMatch = MTetrodes == tetrode;
clusterMatch = MClusters == cluster;


ind = find(dateMatch & tetrodeMatch & clusterMatch);

% Warn if multiple matches. This is a danger because we allow partial date
% matches, and the user may give a degenerate string for the date.
if length(ind) > 1
  warning('findInM:multipleMatches', ...
    'Found %d matching members of M for date %s, tetrode %d, cluster %d', ...
    length(ind), theDate, tetrode, cluster);
end
