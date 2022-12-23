function sdu = twoP_computeZScoreFromSpks(animal,session,overwrite)
% This function is work-in-progress
% Default z-score window: 5 minutes
% Default filter window
% Need to make filter type and filter length input options
%%
% animal = 'Fez59'; session  = '20200608'; overwrite = 1;

if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = 0;
end  

S = twoP_settings;
nMinutes = 5;
filterWindow = 10; filterMethod = 'gauss';
t_int = round(nMinutes*60*S.sRate);

load_filename = 'spks.npy';
load_path = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,load_filename);

save_filename = ['zscore_' num2str(nMinutes) '_minute_filter.mat'];
save_path = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,save_filename);

if isfile(save_path) && (overwrite == 0)
    disp([save_path ' exists, will not replace but will load existing.']);
    sdu = load(save_path,'sdu'); sdu = sdu.sdu;
    return
else
    is = readNPY(load_path); % load inferred spikes
    is = smoothCol(is,2,filterWindow,filterMethod);
    % sdu = (is-movmean(is,t_int,2)) ./ movstd(is,t_int,0,2);
    sdu = bsxfun(@rdivide,is-movmean(is,t_int,2),movstd(is,t_int,0,2));
    sdu = single(sdu);
    save(save_path,'sdu', '-v7.3');
    disp(['ZScore is computed for ' animal ' ' session ' and saved as ' save_path]);
end
