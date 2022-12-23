function twoP_convertSpks2ZScore(animal,session)

S = twoP_settings;
nMinutes = 5;
t_int = round(nMinutes*60*S.sRate);
% file_path = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy\Plex51\imaging\200328\suite2p\plane0\spks.npy';
load_filename = 'spks.npy';
load_path = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,load_filename);
is = readNPY(load_path);

sdu = (is-movmean(is,t_int,2)) ./ movstd(is,t_int,0,2);

save_filename = ['zscore_' num2str(nMinutes) '_minute_filter.mat'];
save_path = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,save_filename);

save(save_path,'sdu');

