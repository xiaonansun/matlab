%rateDisc optogenetics figures
% cPath = '\\grid-hs\churchland_nlsas_data\\data\Behavior_Simon\';
cPath = 'G:\Google Drive\Behavior_Simon\';
groupnames = {'EMX' 'FezF' 'Plexin' 'CSP'};

%% detection performance
% load behavioral data for detection paradigm with optogenetic inactivation
for x = 1 : length(groupnames)
    bhv{x} = rateDisc_loadDetectionBhv(groupnames{x}, cPath, false);
end

% show general effects
for x = 1 : length(groupnames)
    rateDisc_optoBasicFigure(bhv{x},groupnames{x});
end

% control mice
ctrlBhv = rateDisc_loadDetectionBhv('Control', cPath, false);
rateDisc_optoBasicFigure(ctrlBhv, 'Control');

% early vs late in training
rateDisc_optoDiffFigure(bhv,groupnames);

% effects of unilateral stimulation
rateDisc_optoUnilateral(bhv,groupnames);


%% discrimination performance
% load behavioral data for discrimination paradigm with optogenetic inactivation
for x = 1 : length(groupnames)
    fprintf('Loading discrimination data for %s mice\n', groupnames{x});
    bhv{x} = rateDisc_loadDiscriminationBhv(groupnames{x}, cPath, false);
end

% look at discrimination curves
rateDisc_optoDiscCurvesFigure(bhv,groupnames,7);

% reverse correlation
rateDisc_reverseCorrelation(bhv,groupnames,0.125,1);


