

% David Raposo
% Created Mar 2013
% Updated Sep 26, 2015
%
% Added head tracking option. This was previously in a separate function.
% To include head tracking data in the merged data structure you need to
% set the "save_video_data" flag to 1. This is an optional argument and its
% default value is 0. For this to work there must exist a video tracking
% file named "VT1.nvt" inside the raw neural data folder.
% 
% To compute the waveforms of each single unit you must set the
% "save_waveforms" flag to 1. This is an optional argument and its default
% value is 0.
%
% Example calls:
%    all_data = stitchData('2011-07-29_10-31-45','ac004');
%    all_data = stitchData('2011-07-29_10-31-45','ac004', [1:8]);
%    all_data = stitchData('2011-07-29_10-31-45','ac004',[1 3 8]);
% 
% Example call with video tracking data:
%    all_data = stitchData('2012-09-27_16-33-42','ac026', [1:8], 1);
%
% Example call with video tracking data and waveforms of single units:
%    all_data = stitchData('2012-09-27_16-33-42','ac026', [1:8], 1, 1);
%

function [all_data] = stitchData (dirname, ratname, which_tetrodes, save_video_data, save_waveforms)

    if ~exist('save_waveforms', 'var')
        save_waveforms = 0;
    end
    
    if ~exist('save_video_data', 'var')
        save_video_data = 0;
    end
    
    % load behavioral data
    this_date = datestr(dirname(1:10));
   
    % check if this is a mac
    if exist('/Volumes', 'dir')
        % must be a mac
        behav_folder = ['/Volumes/churchland/data/' ratname '/behavior/'];
        dirpath = ['/Volumes/churchland/data/' ratname '/rawneural/' dirname '/'];
    else
        % not a mac, must be a windows computer then
        behav_folder = ['Z:\data\' ratname '\behavior\'];
        dirpath = ['Z:\data\' ratname '\rawneural\' dirname '\'];
    end
    
    if ~exist(behav_folder, 'dir')
        keyboard;
        disp(['Could not find behavior folder: ' behav_folder]);
        exit;
    end
    if ~exist(dirpath, 'dir')
        keyboard;
        disp(['Could not find rawneural folder: ' dirpath]);
        exit;
    end
        
    
    % behav_folder = ['~/Dropbox/BehavData/rats data/' this_date '/'];
    lookfor = [ratname '_' this_date];
    list_files = dir(behav_folder);

    behav_file = [];
    for i = 3:length(list_files)
        if length(list_files(i).name) > length(lookfor)
            if strmatch(list_files(i).name(end-3:end), '.mat')
                if strmatch(list_files(i).name(1:length(lookfor)), lookfor)
                    behav_file = list_files(i).name;
                    break;
                end
            end
        end
    end
    
    if isempty(behav_file)
        disp(['Cannot find behavior file in ' behav_folder]);
        all_data = [];
        return;
    end
    behav_data_temp = load([behav_folder behav_file]);
    behav_data = behav_data_temp.all_data;    

    useTheseFilename = [dirpath 'use_these.mat'];
    if exist(useTheseFilename, 'file')
        disp('Found use_these file');
        use_these = load(useTheseFilename);
        use_these = use_these.use_these;
    end
    
    % if finds both files, will use this last one.
    mclustInfoFilename = [dirpath 'cluster_ratings.mat'];
    if exist(mclustInfoFilename, 'file')
        disp('Found mclust info file');
        mclustInfo = load(mclustInfoFilename);
        use_these = mclustInfo.ratMat;
    end
    
    use_these = sortrows(use_these, [1 2]);
    
    if ~exist('which_tetrodes', 'var')
        which_tetrodes = unique(use_these(:,1));
    end

    events_file = [dirpath 'Events.nev'];
    if ~exist(events_file, 'file')
        disp('Could not find events file.');
        return;
    end
    disp('Reading Events.nev file');
    events = getRawTTLs(events_file);

    

    trials_indexes = find(events(:,2) > 1);
    % This should fix the problem -- David, 12-23-2014
    trials_indexes = trials_indexes(trials_indexes + 1 <= length(events));
    
    valid_check = events(trials_indexes + 1,2) == 1;
    valid_trials_indexes = trials_indexes(valid_check)+1;
    valid_trials_start_time = events(valid_trials_indexes,1);
    valid_trials = events(valid_trials_indexes-1,2);
    
    % check if there are duplicated trial ids
    % cut the repeated ones in the beginning of the session.
    start_here = find(diff(valid_trials) < 0, 1, 'last') + 1;
    if ~isempty(start_here)
        valid_trials = valid_trials(start_here:end);
        valid_trials_start_time = valid_trials_start_time(start_here:end);
    end
        
  
    % check if there is behavioral data for all neuralynx trials
    % this happens for at least one session, not sure why.
    % cut the last trials on neuralynx side.
    keep_these = valid_trials <= length(behav_data);
    valid_trials = valid_trials(keep_these);
    valid_trials_start_time = valid_trials_start_time(keep_these);
    
    
    if exist([dirpath 'VT1.nvt'],'file') < 2
        save_video_data = 0;
        disp('Video file not found');
    end
    
    if save_video_data
        video_file = [dirpath 'VT1.nvt'];
        % the indexes don't really matter, cause we are using mode=1
        % this will extract all the video data and ignore the indexes.
        [video_ts, video_x, video_y, video_angle] = getRawVT(video_file, 23, 240, 1);
        all_video_data = [video_ts' video_x' video_y' video_angle'];
    end
    
  
    
    disp('Merging spikes with behavior:');
    unit_id = 1;

    for i = 1:length(which_tetrodes)
        t = which_tetrodes(i);
        disp(['- Tetrode ' num2str(t)]);
        which_clusters = use_these(use_these(:,1) == t,2);

        ntt_filename = [dirpath 'TT' num2str(t) '.ntt'];

        for j = 1:length(which_clusters)
            c = which_clusters(j);

            % load spikes here
            filename = [dirpath 'TT' num2str(t) '_' num2str(c) '.mat'];
            if exist(filename,'file') == 2
                this_cluster_ts = load(filename);
                spikeTimes = this_cluster_ts.TS;
                
                % The spike times from the .mat file the we are loading are
                % in tenths of miliseconds, but have microsecond precision.
                % The events (start times) are in microseconds already.
                this_iso = use_these(use_these(:,1) == t & use_these(:,2) == c, 3);
               
                these_wave_forms = [];
                if save_waveforms
                    [ntt_spikes ntt_wf] = LoadTT_NeuralynxNT(ntt_filename);
                    % save waveforms
                    [junk,junk2,indexes] = intersect(spikeTimes,ntt_spikes);
                    % [junkTS, wv] = LoadTT_NeuralynxNT(ntt_filename, spikeTimes*100, 1); % Specified spikes
                    % wv: number of spikes * 4 channels * 32 (time points of the waveform)
                    wf = ntt_wf(indexes,:,:);
                    mwf = squeeze(mean(wf,1));
                    stdwf = squeeze(std(wf,1));
                    
                    these_wave_forms.meanWaveform = mwf';
                    these_wave_forms.stdWaveform = stdwf';
                    these_wave_forms.sampleWaveforms = wf(randi(size(wf,1),1,20),:,:);
                end

                behav_data = organize_cluster_trials(behav_data, t, c, unit_id, this_iso, spikeTimes*100, valid_trials, valid_trials_start_time, these_wave_forms);
                unit_id = unit_id + 1;
                
                if save_video_data
                    behav_data = organize_video_trials(behav_data, all_video_data, valid_trials, valid_trials_start_time);
                end

            else
                disp(['Tetrode ' num2str(t) ', Cluster ' num2str(c) ' not found']);
            end
        end
    end    
    
    [behav_data.('clusterInfo')] = deal(use_these);
    
    % clean up the behav_data struct
    behav_data = behav_data(1:end-1);
    if ~isfield(behav_data, 'trialId')
        behav_data = cleanup_behav_data(behav_data); % comment this line if you want the old structure fields
    end
    
    % Fixing missing field
    if ~isfield(behav_data, 'imposedWaitDuration')
        fprintf('Fixing imposedWaitDuration field\n');
        if isfield(behav_data, 'stimWaitDuration')
            [behav_data.imposedWaitDuration] = deal(behav_data(:).stimWaitDuration);
        else
            fprintf('Could not find stimWaitDuration field\n');
        end
    end
    
    for trial = 1:length(behav_data)
        if ~isempty(behav_data(trial).units) && ~isempty(vertcat(behav_data(trial).units(:).spikes))
            behav_data(trial).hasNeuralData = 1;
        else
            behav_data(trial).hasNeuralData = 0;
        end
    end
    
    % add stability info (by Matt)
    % behav_data = markAlldataInstability(behav_data);
    
    % waveforms = all_tetrodes;
    all_data = behav_data;
    % save(['data_20sep12/' ratname '/' ratname '_' datestr(this_date, 'yyyy-mm-dd') '.mat'],'all_data');
    
    disp('Done.');
end


function new_behav_data = organize_cluster_trials (behav_data, tetrode, cluster, unit_id, iso, sTimes, trial_numbers, trial_time_zero, wave_forms)
    sp = 1;
    tr = 1;
    
    pre_trial = -1e6; % in microseconds
    post_trial = 3e6;
    
    nTrials = length(trial_time_zero);
    
    trialStart = trial_time_zero + pre_trial;
    trialEnd = trial_time_zero + post_trial;
    
   
    % Initialize the units field for this unit for each trial
    initUnit.tetrodeNumber = tetrode;
    initUnit.clusterNumber = cluster;
    initUnit.isolation = iso;
    initUnit.spikes = [];
    initUnit.waveforms = wave_forms;
    
    for initTr = 1:length(behav_data)
        behav_data(initTr).units(unit_id) = initUnit;
        behav_data(initTr).preTrial = pre_trial/1000;
        behav_data(initTr).postTrial = post_trial/1000;
        behav_data(initTr).videoTrackingData = [];
        behav_data(initTr).hasVideoTrackingData = 0;
    end

    
    % Check if first spike is a member of the first trial. If not, run forward
    % until we find the first spike.
    while sTimes(sp) < trialStart(tr)
        sp = sp + 1;
    end
    
    while 1
        % Record beginning
        firstSpike = sp;
        
        % Find end of trial
        while sp < length(sTimes) && sTimes(sp) <= trialEnd(tr)
            sp = sp + 1;
        end
        % Record end
        lastSpike = sp - 1;
        
       % save spike times in milisenconds
        behav_data(trial_numbers(tr)).units(unit_id).spikes = (sTimes(firstSpike:lastSpike) - trial_time_zero(tr)) / 1000;
         
        if tr >= nTrials
            break;
        end
        % Search for beginning of next trial
        tr = tr + 1;
        if sTimes(sp) < trialStart(tr)
            % Run forward
            while sp < length(sTimes) && sTimes(sp) < trialStart(tr)
                sp = sp + 1;
            end
            if sp >= length(sTimes)
                break;
            end
        else
            % Run backward
            while sp > 0 && sTimes(sp) >= trialStart(tr)
                sp = sp - 1;
            end
            sp = sp + 1;
        end
    end
    new_behav_data = behav_data;
    
end

function new_behav_data = organize_video_trials (behav_data, video_data, trial_numbers, trial_time_zero)

    frameTimes = video_data(:,1);
    xPosition = video_data(:,2);
    yPosition = video_data(:,3);
    headAngle = video_data(:,4);
    
    vf = 1; % video frame
    tr = 1;
    
    pre_trial = -1e6; % in microseconds
    post_trial = 3e6;
    
    nTrials = length(trial_time_zero);
    
    trialStart = trial_time_zero + pre_trial;
    trialEnd = trial_time_zero + post_trial;
    
    % Check if first spike is a member of the first trial. If not, run forward
    % until we find the first spike.
    while frameTimes(vf) < trialStart(tr)
        vf = vf + 1;
    end
    
    while 1
        % Record beginning
        firstFrame = vf;
        
        % Find end of trial
        while vf <= length(frameTimes) && frameTimes(vf) <= trialEnd(tr)
            vf = vf + 1;
        end
        % Record end
        lastFrame = vf - 1;
        
        % save spike times in milisenconds
        behav_data(trial_numbers(tr)).videoTrackingData.frameTimes = (frameTimes(firstFrame:lastFrame) - trial_time_zero(tr)) / 1000;
        behav_data(trial_numbers(tr)).videoTrackingData.x = xPosition(firstFrame:lastFrame);
        behav_data(trial_numbers(tr)).videoTrackingData.y = yPosition(firstFrame:lastFrame);
        behav_data(trial_numbers(tr)).videoTrackingData.headAngle = headAngle(firstFrame:lastFrame);
        behav_data(trial_numbers(tr)).hasVideoTrackingData = 1;

        
        if tr >= nTrials
            break;
        end
        % Search for beginning of next trial
        tr = tr + 1;
        if frameTimes(vf) < trialStart(tr)
            % Run forward
            while vf <= length(frameTimes) && frameTimes(vf) < trialStart(tr)
                vf = vf + 1;
            end
            if vf >= length(frameTimes)
                break;
            end
        else
            % Run backward
            while vf > 0 && frameTimes(vf) >= trialStart(tr)
                vf = vf - 1;
            end
            vf = vf + 1;
        end
    end
    new_behav_data = behav_data;
end

