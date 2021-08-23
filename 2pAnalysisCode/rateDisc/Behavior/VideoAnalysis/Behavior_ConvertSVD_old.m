function Behavior_ConvertSVD(fPath)
% code to convert behavioral videos from timestamp (.mat) and compressed movie (.mj2) files.

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end
opts.nSVD = 2000;
opts.memLimit = 100; % memory limit for video data in workspace in gigabyte. Use frame averaging to stay within limit.

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

                %% loop through available webcams
                for iCams = 1:nrCams
                    
                    cPath = [fPath aPath  pPath sPath ePath]; %path to current experiment
                    timeFiles = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*frameTimes*_' int2str(iCams) '.mat']); %all files for current cam based on integer before .mat
                    movieFiles = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*video*_' int2str(iCams) '.mj2']); %all files for current cam based on integer before .mat
                    
                    if length(movieFiles) > 5 %don't do conversion if less than 5 movie files are found - probably something was not deleted correctly
                        %% check frameTimes to asses total number of frames in experiment
                        tic
                        fCnt = 0;
                        for iFiles = 1:length({timeFiles.name})
                            cFile = [cPath timeFiles(iFiles).name];
                            load(cFile)
                            fCnt = fCnt + size(frameTimes,1);
                        end
                        
                        %% load single movie and get one frame to compute expected workspace size for whole experiment
                        cFile = [cPath movieFiles(1).name];
                        v = VideoReader(cFile);
                        singleFrame = readFrame(v); clear v
                        info = whos('singleFrame');
                        exptSize = info.bytes * fCnt * 4 / 2^30; %expected size of complete data set in gb
                        frameAvg = ceil(exptSize / opts.memLimit); %average across frames to keep memory usage under control
                        
                        %% load video data and accuumulate in larger array
                        Cnt = 0;
                        mov = zeros(size(singleFrame,1),size(singleFrame,2),floor(fCnt/frameAvg),'single'); %pre-allocate mov array to compute svd later
                        
                        for iFiles = 1:length({movieFiles.name})
                            cFile = [cPath movieFiles(iFiles).name];
                            rawData = squeeze(importdata(cFile));
                            rawData(:,:,1:randi(frameAvg)) = []; %remove a random nr of frames to increase variability in the movie
                            rawData(:,:,floor(size(rawData,3) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above divider
                            rawData = reshape(rawData, size(rawData,1), size(rawData,2), [], frameAvg); %reshape data and keep average in mov array
                            mov(:,:,Cnt + (1:size(rawData,3)),1) = mean(rawData,4);
                            Cnt = Cnt + size(rawData,3);
                            if rem(iFiles,100) == 0
                                fprintf(1, 'Cam%d: Loaded file %d/%d\n',iCams, iFiles,length({timeFiles.name}));
                            end
                        end
                        mov(:,:,Cnt+1:floor(fCnt/frameAvg)) = []; %delete unused frames from mov array
                        
                        %% compute svd
                        tic
                        disp('Computing SVD');
                        opts.nSVD = min(opts.nSVD, size(mov,3));
                        mov       = reshape(mov, [], size(mov,3));
                        COV       = mov' * mov/size(mov,1);
                        totalVar  = sum(diag(COV)); % total variance of data.
                        opts.nSVD = min(size(COV,1)-2, opts.nSVD);
                        [V, Sv]   = svd(COV);
                        
                        V         = V(:, 1:opts.nSVD);
                        Sv        = single(diag(Sv(1:opts.nSVD, 1:opts.nSVD)));
                        U         = single(normc(mov * V));
                        clear COV mov V
                        
                        %% apply to data
                        disp('Compute temporal component V');
                        Cnt = 0;
                        totalFrameTimes = zeros(1,fCnt);
                        V = zeros(opts.nSVD,fCnt,'single');
                        
                        for iFiles = 1:length({movieFiles.name})
                            %get video data and compute temportal component V from spatial component U
                            cFile = [cPath movieFiles(iFiles).name];
                            rawData = squeeze(importdata(cFile));
                            rawData = single(rawData);
                            V(:,Cnt + (1:size(rawData,3))) = U' * reshape(rawData, [], size(rawData,3));
                            delete(cFile);
                            
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
                        U = reshape(U, size(singleFrame,1), size(singleFrame,2), []);
                        disp('Saving data');
                        save([cPath 'SVD_Cam ' int2str(iCams) '.mat'],'V','U', 'totalFrameTimes', 'Sv','totalVar');
                        toc
                    end
                end
            end
        end
    end
end
