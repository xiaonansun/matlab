animal = 'Plex68';

vidDir =  '\\grid-hs\churchland_nlsas_data\BehaviorVideo\';
bhvDir =  '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
imDir =  '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';

vidSessions = dir([vidDir animal filesep 'SpatialDisc\Session Data']);
parfor i = 3:length(vidSessions)
   if  vidSessions(i).isdir ==1
       vidFiles{i} = dir([vidSessions(i).folder filesep vidSessions(i).name]);
   end
end
% vidSessions = dir([vidDir animal filesep 'SpatialDisc\Session Data']);
% sum(cell2mat({vidSessions(3:end).isdir}))