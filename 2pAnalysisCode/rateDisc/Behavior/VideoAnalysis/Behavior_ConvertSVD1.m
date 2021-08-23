function Behavior_ConvertSVD1(fPath)
% code to convert behavioral videos from timestamp (.mat) and compressed movie (.mj2) files.

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end
opts.nSVD = 2000;
opts.memLimit = 60; % memory limit for video data in workspace in gigabyte. Use frame averaging to stay within limit.
opts.varCnt = 100; %number of trials used to compute variance map
opts.eyeFrame = 10; %add some pixels to eyeFrame to ensure pupil is properly captured
opts.eyeThresh = .25;
opts.eyeErrode = 2;
opts.maxFrameCnt = 1000; %max number of frames per trial. Earlier frames will be removed from analysis.

%% folders that contain animal video data
animals = dir(fPath);
animals = animals([animals.isdir] & ~strncmpi('.', {animals.name}, 1));

%% loop through animal data
for iAnimals = 1:length({animals.name})
    
    aPath = [animals(iAnimals).name filesep];
    paradigms = dir([fPath aPath]);
    paradigms = paradigms([paradigms.isdir] & ~strncmpi('.', {paradigms.name}, 1));
    
    %% loop through paradigms
    for iParadigms = 1:length({paradigms.name})
        
        pPath = [paradigms(iParadigms).name filesep];
        sessions = dir([fPath aPath  pPath]);
        sessions = sessions([sessions.isdir] & ~strncmpi('.', {sessions.name}, 1));
        
        %% loop through sessions
        for iSessions = 1:length({sessions.name})
            
            sPath = [sessions(iSessions).name filesep];
            experiments = dir([fPath aPath  pPath sPath]);
            experiments = experiments([experiments.isdir] & ~strncmpi('.', {experiments.name}, 1));
            
            %% loop through experiments
            for iExperiments = 1:length({experiments.name})
                ePath = [experiments(iExperiments).name filesep];
                files = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*frameTimes*.mat']);
                temp = cat(1,files.name);
                ind = strfind(temp(1,:),'_');ind = ind(end);
                nrCams = length(unique(str2num(temp(:,ind+1))));
                disp(['Loading experiment ' experiments(iExperiments).name]);
                
                load([fPath aPath pPath sPath ePath(1:end-1) '.mat']); %load bpod data
                eyePos = SessionData.eyePos; %get eye position and camera for face
                snoutPos = SessionData.snoutPos; %get snout position and camera for face
                
                %% loop through available webcams
                for iCams = 1:nrCams
                    
                    cPath = [fPath aPath  pPath sPath ePath]; %path to current experiment
                    timeFiles = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*frameTimes*_' int2str(iCams) '.mat']); %all files for current cam based on integer before .mat
                    movieFiles = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*video*_' int2str(iCams) '.mj2']); %all files for current cam based on integer before .mat
                    
                    if length(movieFiles) > 5 %don't do conversion if less than 5 movie files are found - probably something was not deleted correctly
                        %% check frameTimes to asses total number of frames in experiment
                        fCnt = 0;
                        timeStamps = cell(1,length({timeFiles.name}));
                        for iFiles = 1:length({timeFiles.name})
                            cFile = [cPath timeFiles(iFiles).name];
                            load(cFile)
                            if size(frameTimes,1) > opts.maxFrameCnt
                                timeStamps{iFiles} = frameTimes(end-opts.maxFrameCnt+1:end);
                            else
                                timeStamps{iFiles} = frameTimes;
                            end
                            fCnt(iFiles) = size(timeStamps{iFiles},1);  %nr of frames per trial
                        end
                        
                        %% load single movie and get one frame to compute expected workspace size for whole experiment
                        if opts.varCnt > length(movieFiles); opts.varCnt = length(movieFiles); end
                        
                        cFile = [cPath movieFiles(1).name];
                        v = VideoReader(cFile);
                        singleFrame = arrayResize(readFrame(v),2); clear v
                        
                        temp = randi(length(movieFiles),opts.varCnt,1);
                        rawVars = zeros([size(singleFrame),opts.varCnt],'single');
                        rawMean = zeros([size(singleFrame),opts.varCnt],'single');
                        for iFiles = 1:opts.varCnt
                            cFile = [cPath movieFiles(temp(iFiles)).name];
                            rawData = squeeze(importdata(cFile));
                            rawVars(:,:,iFiles) = std(single(arrayResize(rawData,2)),[],3);
                            rawMean(:,:,iFiles) = mean(single(arrayResize(rawData,2)),3);
                        end
                        rawVars = mean(rawVars,3); %use variance to select pixels
                        save([cPath 'rawVars_Cam' int2str(iCams) '.mat'],'rawVars');
                        
                        rawMean = mean(rawMean,3); %get mean picture
                        save([cPath 'rawMean_Cam' int2str(iCams) '.mat'],'rawMean');
                        
                        mask = rawVars <= 2; %use variance to select pixels
                        singleFrame = single(singleFrame(~mask)); %get selected pixels in downsampled image
                        info = whos('singleFrame');
                        exptSize = info.bytes * sum(fCnt) / 2^30; %expected size of complete data set in gb
                        frameAvg = ceil(exptSize / opts.memLimit); %average across frames to keep memory usage under control
                        fprintf(1,'MemLimit: %f ; Expected: %f ; FrameAvg: %d \n', opts.memLimit,exptSize,frameAvg);
                        
                        %% load video data and accuumulate in larger array
                        Cnt = 0;
                        mov = zeros(sum(floor(fCnt/frameAvg)),size(singleFrame,1),'single'); %pre-allocate mov array to compute svd later
                        
                        for iFiles = 1:length({movieFiles.name})
                            cFile = [cPath movieFiles(iFiles).name];
                            rawData = squeeze(importdata(cFile));
                            if size(rawData,3) > opts.maxFrameCnt
                                rawData = rawData(:,:,end-opts.maxFrameCnt+1:end);
                            end
                            
                            %% get trace for eye video and compute pupil size
                            if iCams == eyePos(end)
                                eyeTrace = rawData((eyePos(2)-eyePos(3)-(opts.eyeFrame/2)-1) + (1:eyePos(3)*2+opts.eyeFrame+1), (eyePos(1)-eyePos(3)-(opts.eyeFrame/2)-1) + (1:eyePos(3)*2+opts.eyeFrame+1), :);
                                eyeTrace = mat2gray(eyeTrace);
                                
                                eyeRad = zeros(1,size(eyeTrace,3));
                                eyeCom = zeros(size(eyeTrace,3),2);
                                for iFrames = 1:size(eyeTrace,3)
                                    eyeVars = Behavior_EyeCheck(eyeTrace(:,:,iFrames),opts.eyeThresh,opts.eyeErrode); %get pupild data
                                    eyeRad(iFrames) = eyeVars.rad; %pupil radius
                                    eyeCom(iFrames,:) = eyeVars.com; %pupil position
                                end
                                eyeTrace = im2uint8(reshape(eyeTrace,size(eyeTrace,1),size(eyeTrace,2),1,[]));
                                
                                v = VideoWriter([cPath 'eyeTrace_' int2str(iFiles) '.mj2'],'Archival'); %write small eye video
                                open(v);writeVideo(v,eyeTrace);close(v);
                                
                                % compute snout motion
                                snoutMotion = rawData((snoutPos(2)-snoutPos(3)-1) + (1:snoutPos(3)*2+1), (snoutPos(1)-snoutPos(3)-1) + (1:snoutPos(3)*2+1), :);
                                snoutMotion = diff(snoutMotion,1,3);
                                snoutMotion = reshape(snoutMotion,[],size(rawData,3)-1);
                                snoutMotion = mean(snoutMotion);
                                snoutMotion = [snoutMotion(1) snoutMotion];
                                
                                %coompute face motion
                                faceMotion = diff(rawData,1,3);
                                faceMotion = reshape(faceMotion,[],size(rawData,3)-1);
                                faceMotion = mean(faceMotion);
                                faceMotion = [faceMotion(1) faceMotion];
                                
                                cTimes = ((timeStamps{iFiles} - SessionData.TrialStartTime(iFiles)) * 86.4*1e3);
                                save([cPath 'faceVars_' int2str(iFiles) '.mat'],'eyeRad','eyeCom','snoutMotion','faceMotion','cTimes');
                            end
                            
                            %% merge data in mov array for svd compression
                            rawData = single(arrayResize(rawData,2));
                            rawData = arrayShrink(rawData,mask,'merge');
                            rawData(:,1:randi(frameAvg)) = []; %remove a random nr of frames to increase variability in the movie
                            rawData(:,floor(size(rawData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above divider
                            
                            rawData = reshape(rawData, size(rawData,1), [], frameAvg); %reshape data and keep average in mov array
                            mov(Cnt + (1:size(rawData,2)), :) = mean(rawData,3);
                            Cnt = Cnt + size(rawData,2);
                            
                            if rem(iFiles,50) == 0
                                fprintf(1, 'Cam%d: Loaded file %d/%d\n',iCams, iFiles,length({timeFiles.name}));
                            end
                        end
%                         mov(:,Cnt+1:end) = []; %delete unused frames from mov array. This can take a lot of processing time.
                        
                        %% compute svd
                        tic
                        disp('Computing SVD');
                        opts.nSVD = min(opts.nSVD, Cnt);
                        COV       = double(mov(1:Cnt,:)' *  mov(1:Cnt,:)/size(mov,1));
                        totalVar  = sum(diag(COV)); % total variance of data.
                        opts.nSVD = min(size(COV,1)-2, opts.nSVD);
                        
                        [U, Sv]   = eigs(COV, opts.nSVD, 'la');
                        Sv        = single(diag(Sv(1:opts.nSVD, 1:opts.nSVD)));
                        
                        clear COV mov
                        
                        %% apply to data
                        disp('Compute temporal component V');
                        Cnt = 0;
                        totalFrameTimes = zeros(1,sum(fCnt));
                        V = zeros(opts.nSVD,sum(fCnt),'single');
                        
                        for iFiles = 1:length({movieFiles.name})
                            %get video data and compute temportal component V from spatial component U
                            cFile = [cPath movieFiles(iFiles).name];
                            rawData = squeeze(importdata(cFile));
                            if size(rawData,3) > opts.maxFrameCnt
                                rawData = rawData(:,:,end-opts.maxFrameCnt+1:end);
                            end
                            rawData = single(arrayResize(rawData,2));
                            rawData = arrayShrink(rawData,mask,'merge');
                            
                            V(:,Cnt + (1:size(rawData,2))) = U' * rawData;
%                             delete(cFile);
                            
                            % combine all frameTimes into single variable
                            cFile = [cPath timeFiles(iFiles).name];
                            load(cFile)
                            frameTimes = frameTimes; %convert to ms
                            totalFrameTimes(1,Cnt + (1:length(frameTimes))) = frameTimes;
                            
                            Cnt = Cnt + size(rawData,3);
                            
                            if rem(iFiles,100) == 0
                                fprintf(1, 'Cam%d: Loaded file %d/%d\n',iCams, iFiles,length({timeFiles.name}));
                            end
                        end
                        U = reshape(U, size(mask,1), size(mask,2), []);
                        disp('Saving data');
                        save([cPath 'SVD_Cam ' int2str(iCams) '.mat'],'V','U', 'totalFrameTimes', 'Sv','totalVar');
                        clear V U
                        toc
                    end
                end
                cTime = now;
                save([cPath 'SVD_Complete.mat'],'cTime');
            end
        end
    end
end
