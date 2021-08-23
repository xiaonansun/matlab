function out = rateDisc_shrinkBhv(bhv)
% code to reduce bhv array to fields that are actually needed. Otherwise it
% becomes very large and slows down the analysis.

%% fields to keep
keepers = {'optoPower' 'stimLocation' 'date' 'SessionNr' 'AnimalID' 'optoDur' 'optoType' 'optoSide' ...
           'StimType' 'DistStim' 'CorrectSide'};

%% keep some fields. use all logicals.
allFields = fields(bhv);
for iFields = 1 : length(allFields)
    cData = bhv.(allFields{iFields});
    
    if islogical(cData)
        out.(allFields{iFields}) = cData;
        
    elseif any(ismember(keepers, allFields{iFields}))
        out.(allFields{iFields}) = single(cData);
    end
end
