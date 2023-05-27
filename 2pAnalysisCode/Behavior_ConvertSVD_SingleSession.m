function Behavior_ConvertSVD_SingleExp(cPath,newRun,invertPupil)
%%
% newRun = false;
% inverPupil = true;
% cPath = '\\grid-hs\churchland_nlsas_data\BehaviorVideo\Fez61\SpatialDisc\Session Data\Fez61_SpatialDisc_Jun04_2020_Session2';

idx_path = regexp(cPath,filesep);
ePath = cPath(idx_path(end)+1:end);
fPath = cPath(idx_path(1):idx_path(6));
pPath = cPath(idx_path(6)+1:idx_path(7));
sPath = cPath(idx_path(7)+1:idx_path(8));

if ePath(end) ~= filesep
    ePath = [ePath filesep];
end

%% loop through experiments
%         for iExperiments = 1:length({experiments.name})
bhvOpts.nSVD = 500;
bhvOpts.memLimit = 50; % memory limit for video data in workspace in gigabyte. Use frame averaging to stay within limit.
bhvOpts.varCnt = 50; %number of trials used to compute variance map
bhvOpts.eyeFrame = 10; %add some pixels to eyeFrame to ensure pupil is properly captured
bhvOpts.eyeThresh = .25;
bhvOpts.eyeErrode = 2;
bhvOpts.maxFrameCnt = 500; %max number of frames per trial. Earlier frames will be removed from analysis.
bhvOpts.memLimit = 32; % memory limit for video data in workspace in gigabyte. Determines how many segments can be computed at once.
bhvOpts.segSize = 4; %segment size relative to full video. Minimum is 2. Will create segSize^2 segments and create svd for each. Second SVD combines all segments into one V/U.
bhvOpts.targRate = 30; %rate of widefield imaging in Hz. For combined SVD, video data is resampled to match targRate.
bhvOpts.preStimDur = 3; %baseline duration for trial-aligned video data
bhvOpts.postStimDur = 115/bhvOpts.targRate; % Duration of trial after stimulus onset in seconds
bhvOpts.framesPerTrial = round((bhvOpts.preStimDur + bhvOpts.postStimDur) * bhvOpts.targRate); %amount of frames per trial in trial-aligned video data
bhvOpts.trialDur = bhvOpts.framesPerTrial / bhvOpts.targRate; %amount of frames per trial in trial-aligned video data
bhvOpts.smth = 5; % filter length when smoothing video data before combined svd
movType = '.mj2';

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

% ePath = [experiments(iExperiments).name filesep];
files = dir([fPath pPath sPath ePath '*' ePath(1:end-1) '*frameTimes*.mat']);
if ~isempty(files)
    temp = cat(1,files.name);
    ind = strfind(temp(1,:),'_');ind = ind(end);
    nrCams = length(unique(str2num(temp(:,ind+1))));
end
svdCheck = dir([fPath pPath sPath ePath '*SVD_CombinedSegments.mat']);
if isempty(svdCheck)
    svdCheck = dir([fPath pPath sPath ePath '*SVD_Complete.mat']);
end

%%
if (length(svdCheck) == 2) && ~newRun
    disp(['Found SVDcomplete in experiment: ' cPath ' - skipped']);
elseif isempty(files)
    disp(['No video data found in experiment: ' cPath ' - skipped']);
else
    try
        %%
        disp(['Loading experiment: ' cPath]);
        load([fPath pPath sPath ePath(1:end-1) '.mat']); %load bpod data
        eyePos = []; snoutPos = []; nosePos = [];
        
        if isfield(SessionData,'eyePos')
            eyePos = SessionData.eyePos; %get eye position and camera for face
        end
        
        if isfield(SessionData,'snoutPos')
            snoutPos = SessionData.snoutPos; %get snout position and camera for face
        end
        
        if isfield(SessionData,'nosePos')
            nosePos = SessionData.nosePos; %get snout position and camera for face
        end
        
        trTime=cell(nrCams,1); % generate data to keep track of trials numbers for timestamp
        trMov=cell(nrCams,1); % generate data to keep track of trials numbers for movie data
        
        %% loop through available webcams
        for iCams = 1:nrCams
            %%
            cPath = [fPath pPath sPath ePath]; %path to current experiment
            timeFiles = dir([fPath pPath sPath ePath '*' ePath(1:end-1) '*frameTimes*_' int2str(iCams) '.mat']); %all files for current cam based on integer before .mat
            movieFiles = dir([fPath pPath sPath ePath '*' ePath(1:end-1) '*Video*_' int2str(iCams) movType]); %all files for current cam based on integer before .mat
            trials = 1:length(movieFiles);
            
            if length(movieFiles) > 1 %don't do conversion if less than 25 movie files are found - probably something was not deleted correctly
                
                %% check if all files have trialcount that matches Sessiondata. Very rarely there can be one too many.
                for iFiles = 1 : length(timeFiles)
                    cTrial = textscan(timeFiles(iFiles).name, '%s', 'Delimiter' , '_');
                    cTrial = str2double(cTrial{1}{end-1});
                    trTime{iCams}(iFiles) = cTrial;
                    if cTrial > length(SessionData.Rewarded)
                        timeFiles(iFiles) = [];
                        warning('Found timestamp file with trial higher trialcount as contained in SessionData. Rejected!')
                    end
                end
                
                for iFiles = 1 : length(movieFiles)
                    cTrial = textscan(movieFiles(iFiles).name, '%s', 'Delimiter' , '_');
                    cTrial = str2double(cTrial{1}{end-1});
                    trMov{iCams}(iFiles) = cTrial;
                    if cTrial > length(SessionData.Rewarded)
                        movieFiles(iFiles) = [];
                        warning('Found movie file with trial higher trialcount as contained in SessionData. Rejected!')
                    end
                end
                timeTemp = false(1,max(trTime{iCams})); timeTemp(trTime{iCams}) = true; trTime{iCams} = timeTemp; clear timeTemp;
                movTemp = false(1,max(trMov{iCams})); movTemp(trMov{iCams}) = true; trMov{iCams} = movTemp; clear movTemp;
                useTrial(iCams,:) = trTime{iCams} & trMov{iCams};
                
                %% check frameTimes and save totalFrameTimes for each camera
                frameCnt = 0;
                timeStamps = cell(1,length({timeFiles.name}));
                for iFiles = 1:length(timeFiles)
                    load([cPath timeFiles(iFiles).name]) % load frame times
                    
                    % check if timestamps are shifted by an hour. Apparently that can happen.
                    timeCheck = (SessionData.TrialStartTime(iFiles)*86400) - (frameTimes(1:10) * 86400); %time difference between first acquired frame and onset of current trial
                    if any(timeCheck > 3540 & timeCheck < 3660) %timeshift by one hour (+- 10seconds)
                        if iFiles == 1
                            warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
                        end
                        frameTimes = frameTimes + (1/24); %add one hour
                    elseif any(timeCheck > 30)
                        error(['Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds in file: ' timeFiles(iFiles).name])
                    end
                    
                    if size(frameTimes,1) > bhvOpts.maxFrameCnt
                        timeStamps{iFiles} = frameTimes(end-bhvOpts.maxFrameCnt+1:end);
                    else
                        timeStamps{iFiles} = frameTimes;
                    end
                    frameCnt(iFiles) = size(timeStamps{iFiles},1);  %nr of frames per trial
                end
                totalFrameTimes = cat(1,timeStamps{:});
                save([cPath 'SVD_Cam' int2str(iCams) '-frameTimes.mat'],'totalFrameTimes','frameCnt');
                
                %% get single frame to assess video size
                cFile = [cPath movieFiles(1).name];
                v = VideoReader(cFile);
                singleFrame = single(arrayResize(readFrame(v),2)); clear v
                if size(singleFrame,3) > 1; singleFrame = singleFrame(:,:,1); end
                
                blockSize = ceil(size(singleFrame)/bhvOpts.segSize);
                ind = im2col(reshape(1:numel(singleFrame),size(singleFrame)),blockSize,'distinct'); %build index for diffent image blocks
                save([cPath 'segInd' int2str(iCams) '.mat'],'ind'); %save index for later reconstruction
                
                %% Check all videos first to ensure none are corrupt
                %                 corruptVideos = false(nrCams,length(movieFiles));
                %                 corruptErrors = cell(nrCams,length(movieFiles));
                %                 movieFileNameList = {movieFiles.name};
                %                 movieFolder = movieFiles(1).folder;
                %                 parfor y = 1:length(movieFiles)
                %                     cFile = fullfile(movieFolder,movieFileNameList{y});
                %                     v = VideoReader(cFile);
                %                     disp(cFile)
                %                     try
                %                     while hasFrame(v)
                %                         frame = readFrame(v);
                %                     end
                %                     catch ME
                %                         corruptVideos(iCams,y) = true;
                %                         corruptErrors{iCams,y} = ME;
                %                     end
                % %                     clear v;
                %                 end
                
                %% determine how often to load raw data given available memory
                info = whos('singleFrame');
                exptSize = (info.bytes * sum(frameCnt) / 2^30) / size(ind,2); %expected size in gb for each video segment
                segPerRun = floor(bhvOpts.memLimit / exptSize); %number of segments per run
                segRuns = ceil(size(ind,2) / segPerRun); %required amount of loading raw data to stay within memory constraints
                fprintf(1, 'Cam%d: Expected datasize: %fgb. Reload data %d times.\n',iCams, exptSize, segRuns);
                segCnt = 0;
                
                %% compute svd for each video segment
                for iSegs = 1 : segRuns
                    if iSegs == segRuns %for last run, check if segPerRun should be reduced
                        segPerRun = size(ind,2) - segCnt;
                    end
                    for x = 1 : segPerRun
                        mov{x} = zeros(sum(frameCnt),sum(ind(:,segCnt + x) <= numel(singleFrame)),'single'); %pre-allocate mov array to compute svd later. Allows for different segment sizes if index is larger as a single frame.
                    end
                    
                    Cnt = 0;
                    for iFiles = 1:length(movieFiles)
                        %% 
                        % iFiles = 79;
                        try
                            errorList{iCams,iFiles} = []; % track error messages of each trial
                            skipTrials(iCams,iFiles) = false; % track indices of corrupt trial videls
                            cFile = [cPath movieFiles(iFiles).name];
                            rawData = squeeze(importdata(cFile));
                            if size(rawData,3) == 3; rawData = squeeze(rawData(:,:,1,:));end
                            
                            if size(rawData,3) > bhvOpts.maxFrameCnt
                                rawData = rawData(:,:,end-bhvOpts.maxFrameCnt+1:end);
                            end
                            
                            %% get trace for eye video and compute pupil size
                            if iSegs == 1
                                
                                cTimes = ((timeStamps{iFiles} - SessionData.TrialStartTime(iFiles)) * 86400); %camera time stamps, relative to trial onset
                                if iCams == eyePos(end) %face camera. Isolate movie for the eye and save snout, nose and facemotion
                                    eyeTrace = rawData((eyePos(2)-eyePos(3)-(bhvOpts.eyeFrame/2)-1) + (1:eyePos(3)*2+bhvOpts.eyeFrame+1), (eyePos(1)-eyePos(3)-(bhvOpts.eyeFrame/2)-1) + (1:eyePos(3)*2+bhvOpts.eyeFrame+1), :);
                                    eyeTrace = mat2gray(eyeTrace);
                                    
                                    eyeTrace = im2uint8(reshape(eyeTrace,size(eyeTrace,1),size(eyeTrace,2),1,[]));
                                    v = VideoWriter([cPath 'eyeTrace_' int2str(iFiles) '.mj2'],'Archival'); %write small eye video
                                    open(v); writeVideo(v,eyeTrace); close(v);
                                    
                                    % compute snout motion
                                    idx1 = (snoutPos(2)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
                                    idx2 = (snoutPos(1)-snoutPos(3)-1) + (1:snoutPos(3)*2+1);
                                    if any(idx1 > size(rawData,1)) || any(idx2 > size(rawData,2))
                                        if iFiles == 1
                                            warning('Snout index is too large. Clipped to fit into image.');
                                        end
                                        idx1(idx1 > size(rawData,1)) = [];
                                        idx2(idx2 > size(rawData,2)) = [];
                                    end
                                    snoutMotion = rawData(idx1, idx2, :);
                                    snoutMotion = diff(snoutMotion,1,3);
                                    snoutMotion = reshape(snoutMotion,[],size(rawData,3)-1);
                                    snoutMotion = mean(snoutMotion);
                                    snoutMotion = [snoutMotion(1) snoutMotion];
                                    
                                    %coompute face motion
                                    faceMotion = diff(rawData,1,3);
                                    faceMotion = reshape(faceMotion,[],size(rawData,3)-1);
                                    faceMotion = mean(faceMotion);
                                    faceMotion = [faceMotion(1) faceMotion];
                                    
                                    %compute nose motion
                                    idx1 = (nosePos(2)-nosePos(3)-1) + (1:nosePos(3)*2+1);
                                    idx2 = (nosePos(1)-nosePos(3)-1) + (1:nosePos(3)*2+1);
                                    if any(idx1 > size(rawData,1)) || any(idx2 > size(rawData,2))
                                        if iFiles == 1
                                            warning('Nose index is too large. Clipped to fit into image.');
                                        end
                                        idx1(idx1 > size(rawData,1)) = [];
                                        idx2(idx2 > size(rawData,2)) = [];
                                    end
                                    noseMotion = rawData(idx1, idx2, :);
                                    noseMotion = diff(noseMotion,1,3);
                                    noseMotion = reshape(noseMotion,[],size(rawData,3)-1);
                                    noseMotion = mean(noseMotion);
                                    noseMotion = [noseMotion(1) noseMotion];
                                    
                                    eyeVars = []; %pupil is now computed later
                                    save([cPath 'faceVars_' int2str(iFiles) '.mat'],'eyeVars','snoutMotion','faceMotion','noseMotion','cTimes');
                                    
                                else %body camera, get body motion and save
                                    bodyMotion = diff(rawData,1,3);
                                    bodyMotion = reshape(bodyMotion,[],size(rawData,3)-1);
                                    bodyMotion = mean(bodyMotion);
                                    bodyMotion = [bodyMotion(1) bodyMotion];
                                    save([cPath 'bodyVars_' int2str(iFiles) '.mat'],'bodyMotion','cTimes');
                                end
                            end
                            
                            %% merge data in mov array for svd compression
                            rawData = arrayResize(rawData,2);
                            rawData = reshape(rawData, [], size(rawData,3)); %reshape data
                            
                            for x = 1 : segPerRun
                                cInd = ind(:,segCnt + x);
                                cInd(cInd > size(rawData,1)) = []; %ensure index does not go beyond available pixels
                                mov{x}(Cnt + (1:size(rawData,2)), :) =  rawData(cInd,:)';  %get current segment
                            end
                            Cnt = Cnt + size(rawData,2);
                            
                            if rem(iFiles,50) == 0
                                fprintf(1, 'Cam%d: Loaded file %d/%d\n',iCams, iFiles,length(trials));
                            end
                        catch ME
                            errorList{iCams,iFiles} = ME;
                            skipTrials(iCams,iFiles) = true;
                        end
                    end
                    
                    if Cnt ~= size(mov{1},1)
%                         error(['Something is wrong with mov array. Cnt = ' num2str(Cnt) '; size(mov{1}) = ' num2str(size(mov{1}))])
                        warning(['Something is wrong with mov array. Cnt = ' num2str(Cnt) '; size(mov{1}) = ' num2str(size(mov{1})) '. Some video files may be corrupt.'])
                    end
                    
                    %% compute svd for current segments
                    for x = 1 : segPerRun
                        tic; fprintf(1, 'Computing SVD - Cam %d, Segment %d\n',iCams,segCnt + x);
                        [V,S,U]  = fsvd(mov{x}, bhvOpts.nSVD);
                        V = V * S;
                        save([cPath 'SVD_Cam' int2str(iCams) '-Seg' int2str(segCnt + x) '.mat'],'V','U');
                        toc
                    end
                    clear mov V U frameTimes
                    segCnt = segCnt + segPerRun;
                end
            end
        end
        
        %% load svd data and merge different segments together
        idx_good_trials = prod(~skipTrials,1); %% RS: index of non-corrupt trials
        
        allV = cell(1,nrCams);
        allTimes = cell(1,nrCams);
        for iCams = 1:nrCams %get segments and frameTimes for each cam
            load([cPath 'SVD_Cam' int2str(iCams) '-frameTimes.mat'],'totalFrameTimes','frameCnt');
            
            allTimes{iCams} = totalFrameTimes; %get timestamps for current cam
            allFrameCnt{iCams} = frameCnt;
            bhvOpts.frameRate(iCams) = round(1/median(diff(allTimes{iCams} * 86400))); %determine frameRate for current cam
            
            % Get all available save segment indices % added by RS
            segFiles = dir([cPath filesep 'SVD_Cam' int2str(iCams) '-Seg*.mat']); %all files for current cam based on integer before .mat

            % load compressed video segments
%             allV{iCams} = cell(1,size(ind,2)); % commented by RS
            allV{iCams} = cell(1,length(segFiles));
%             for iSegs = 1 : size(ind,2)
            for iSegs = 1 : length(segFiles)
                load_file_name = ['SVD_Cam' int2str(iCams) '-Seg' int2str(iSegs) '.mat'];
                load_file_path = [cPath load_file_name];
                disp(['Loading ' load_file_name '...'])
                data = load(load_file_path,'V');
%                 data = load(fullfile(segFiles(iSegs).folder,segFiles(iSegs).name),'V');
                allV{iCams}{iSegs} = data.V; clear data load_file_path
                disp([load_file_name ' added.'])
            end
            allV{iCams} = cat(2,allV{iCams}{:}); % combine in on larger array
        end
        %%
        idx_trial_skip = find(~idx_good_trials);
        for iCams = 1:nrCams
            for iSkip = 1:length(idx_trial_skip)
                skipFrames{iCams,iSkip} = sum(allFrameCnt{iCams}(1:idx_trial_skip(iSkip)-1))+1:sum(allFrameCnt{iCams}(1:idx_trial_skip(iSkip)));
            end
        end
        %% go through each trial and build combined dataset for each cam
        for x = 1 : nrCams
            alignV{x} = cell(1,size(SessionData.RawEvents.Trial,2));
        end
        
        for iTrials = 1:size(SessionData.RawEvents.Trial,2)
            try
                stimTime(iTrials) = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
            catch
                stimTime(iTrials) = NaN;
            end
            
            for iCams = 1:nrCams
                try
                    trialOn = (SessionData.TrialStartTime(iTrials) * 86400) + (stimTime(iTrials) - bhvOpts.preStimDur); % trial onset times
                    cTimes = allTimes{iCams} * 86400 - trialOn; % get frame times for current cam
                    trialOn = find(cTimes > 0,1); %trial on in frames
                    trialOff = find(cTimes > bhvOpts.trialDur,1); %trial off in frames
                    
                    trialTime = cTimes(trialOn:trialOff) - cTimes(trialOn);
                    idx = trialTime < bhvOpts.trialDur; %don't use late frames
                    trialTime = trialTime(idx);
                    timeLeft = bhvOpts.trialDur - trialTime(end); %check if there is missing time at the end of a trial
                    
                    cData = allV{iCams}(trialOn:trialOff,:); %get data based on time vector matching duration of a trial
                    cData = cData(idx, :); %don't use late frames
                    
                    if (timeLeft < bhvOpts.trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
                        addTime = trialTime(end) + (1/bhvOpts.targRate : 1/bhvOpts.targRate : timeLeft + 1/bhvOpts.targRate); %add some dummy times to make complete trial
                        trialTime = [trialTime' addTime];
                        cData = [cData; repmat(cData(1,:),length(addTime),1)];
                    end
                    
                    if any(diff(trialTime) > 0.5) %dont use current trial if frames are no more than 0.5s apart
                        error;
                    end
                    
                    %create time vector to capture potential irregular sampling
                    timePad = 1/bhvOpts.targRate:1/bhvOpts.targRate:1/bhvOpts.targRate*21;
                    trialTime = trialTime'+timePad(end)+1/bhvOpts.frameRate(iCams);
                    trialTime = [timePad trialTime' timePad+trialTime(end)];
                    cData = [repmat(cData(1,:),21,1); cData; repmat(cData(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
                    cData = resample(double(cData), trialTime, bhvOpts.targRate); %resample to set target framerate
                    cData = conv2(cData',ones(1,bhvOpts.smth)/bhvOpts.smth,'same')'; %smooth trace with moving average of 'smth' points
                    alignV{iCams}{iTrials} = cData(22 : 21 + bhvOpts.framesPerTrial,:); %remove padds and select trial-relevant frames
                    
                catch
                    alignV{iCams}{iTrials} = NaN(bhvOpts.framesPerTrial,size(allV{iCams},2),'single'); %dont use this trial
                end
            end
        end
        clear allV
        
        % combined V segments from all cams
        for x = 1 : nrCams
            alignV{x} = cat(1,alignV{x}{:}); %combined V matrix
        end
        alignV = cat(2,alignV{:}); %combined V matrix
        camIdx = ~isnan(mean(alignV,2)); %index for trials where bhv frames were available
        
        % create second SVD
        tic; disp('Compute combined SVD.');
        [V,S,vidU]  = fsvd(alignV(camIdx,:), bhvOpts.nSVD); clear alignV
        V = V * S;
        
        %zscore video regressor, save mean and std for later reconstruction
        meanV = mean(V);
        stdV = std(V);
        V = bsxfun(@minus,V,meanV);
        V = bsxfun(@rdivide,V,stdV);
        
        % rebuild
        vidV = NaN(length(camIdx),bhvOpts.nSVD,'single');
        vidV(camIdx,:) = V; clear V
        save([cPath 'SVD_CombinedSegments.mat'], 'vidV', 'vidU', 'meanV' ,'stdV','allTimes');
        bhvOpts.fPath = cPath;
        save([cPath 'bhvOpts.mat'], 'bhvOpts');
        disp('Done.'); toc;
        
        % compute filtered pupil data
        Behavior_computePupilVars(cPath,true,invertPupil);
        
        % compute svd for absolute motion
        Behavior_AbsoluteMotionSVD([],[],[],bhvOpts,movType);
        
    catch ME
        disp('!!Conversion failed!!');
        disp(ME.message);
    end
end