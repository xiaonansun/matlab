function Behavior_computeLicktrace(fPath)
% code to convert behavioral videos from timestamp (.mat) and compressed movie (.mj2) files.

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end
opts.nSVD = 200;
opts.memLimit = 60; % memory limit for video data in workspace in gigabyte. Use frame averaging to stay within limit.
opts.varCnt = 50; %number of trials used to compute variance map
opts.eyeThresh = .25;
opts.eyeErrode = 2;
opts.maxFrameCnt = 1000; %max number of frames per trial. Earlier frames will be removed from analysis.
Paradigm = 'SpatialDisc';
botCam = 2;
mouthFrame = 50;

%% loop through paradigms
sessions = dir(fPath);
sessions = sessions([sessions.isdir] & ~strncmpi('.', {sessions.name}, 1));

%% loop through sessions
for iSessions = 1:length({sessions.name})
    % for iSessions = 2
    tic
    sPath = [sessions(iSessions).name filesep];
    files = dir([fPath sPath '*' Paradigm '*.mat']);
    
    if ~isempty(files)
        load([fPath sPath files(1).name]); %load bpod data
        movieFiles = dir([fPath sPath 'BehaviorVideo' filesep '*Video*_2.mj2']);
        mouthFiles = dir([fPath sPath 'BehaviorVideo' filesep '*mouthTrace*.mj2']);
        timeFiles = dir([fPath sPath 'BehaviorVideo' filesep '*frameTimes*_2.mat']); %time stamps for bottom cam
        
        if ~isempty(movieFiles)
            if ~isfield(SessionData,'mouthPos')
                cFile = [fPath sPath 'BehaviorVideo' filesep movieFiles(1).name];
                v = VideoReader(cFile);
                singleFrame = readFrame(v);
                
                clear v
                h = figure('units','normalized','outerposition',[0 0 1 1]);
                imagesc(singleFrame); colormap gray; axis image
                mouthPos = ginput(1);
                mouthPos = [round(mouthPos) 50 botCam];
                close(h);
                SessionData.mouthPos = mouthPos;
                save([fPath sPath files(1).name],'SessionData'); %save bpod data
            else
                mouthPos = SessionData.mouthPos;
            end
            mouthInd = [(mouthPos(2)-mouthFrame-1) + (1:mouthFrame*2+1);(mouthPos(1)-mouthFrame-1) + (1:mouthFrame*2+1)];
            
            lickFramesL = cell(1,length(movieFiles));
            lickFramesR = cell(1,length(movieFiles));
            stimFrames = cell(1,length(movieFiles));
            diffStimFrames = cell(1,length(movieFiles));
            frameCnt = zeros(1,length(movieFiles));
            
            %% get trace for eye video and compute pupil size
            for iFiles = 1:length(movieFiles)
                
                if exist([fPath sPath 'BehaviorVideo' filesep 'mouthTrace_' int2str(iFiles) '.mj2'], 'file') ~= 2
                    %                 if true
                    cFile = [fPath sPath 'BehaviorVideo' filesep movieFiles(iFiles).name];
                    rawData = squeeze(importdata(cFile));
                    if size(rawData,3) > opts.maxFrameCnt
                        rawData = rawData(:,:,end-opts.maxFrameCnt+1:end);
                    end
                    mouthTrace = rawData((mouthPos(2)-mouthFrame-1) + (1:mouthFrame*2+1), (mouthPos(1)-mouthFrame-1) + (1:mouthFrame*2+1), :);
                    mouthTrace = mat2gray(mouthTrace);
                    
                    mouthTrace = im2uint8(reshape(mouthTrace,size(mouthTrace,1),size(mouthTrace,2),1,[]));
                    v = VideoWriter([fPath sPath 'BehaviorVideo' filesep 'mouthTrace_' int2str(iFiles) '.mj2'],'Archival'); %write small mouth video
                    open(v);writeVideo(v,mouthTrace);close(v);
                    mouthTrace = squeeze(mouthTrace);
                    
                else
                    mouthTrace = squeeze(importdata([fPath sPath 'BehaviorVideo' filesep 'mouthTrace_' int2str(iFiles) '.mj2']));
                end
                
                % get timestamps for bottom camera
                load([fPath sPath 'BehaviorVideo' filesep timeFiles(iFiles).name]);
                if size(frameTimes,1) > opts.maxFrameCnt
                    frameTimes = frameTimes(end-opts.maxFrameCnt+1:end);
                end
                botTimes = (frameTimes - SessionData.TrialStartTime(iFiles)) * 86.4*1e3 ;
                frameCnt(iFiles) = length(botTimes);
                
                % check for stimulus period
                stimOn = SessionData.RawEvents.Trial{iFiles}.States.PlayStimulus(1);
                stimOff = SessionData.RawEvents.Trial{iFiles}.States.MoveSpout(1); %onset of decision period
                stimTime = botTimes > stimOn & botTimes < stimOff; %index for stimulus period
                stimFrames{iFiles} = mouthTrace(:,:,stimTime);
                
                %subtract average
                mouthTrace = bsxfun(@minus,mouthTrace,uint8(mean(mouthTrace,3)));
                
                % apply temporal filter to data to enhance tongue frames
                a = cat(3,zeros(size(mouthTrace,1),size(mouthTrace,2)), diff(mouthTrace,1,3));
                b = cat(3,flip(uint8(diff(flip(mouthTrace,3),1,3)),3),zeros(size(mouthTrace,1),size(mouthTrace,2)));
                mouthTrace = mean(cat(4,a,b),4);
                diffStimFrames{iFiles} = mouthTrace(:,:,stimTime);

                % check for baseline period
                baseOn = SessionData.RawEvents.Trial{iFiles}.States.WaitBeforeLever(end)-1;
                baseOff = SessionData.RawEvents.Trial{iFiles}.States.WaitBeforeLever(end);
                baseTime = botTimes > baseOn & botTimes < baseOff; %index for baseline period
                baseFrames{iFiles} = mouthTrace(:,:,baseTime);

                % find lick frames
                lLickOn = [];
                if isfield(SessionData.RawEvents.Trial{iFiles}.Events,'Port1In') %check for right licks
                    licks = SessionData.RawEvents.Trial{iFiles}.Events.Port1In;
                    licks(licks < SessionData.RawEvents.Trial{iFiles}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
                    
                    for iLicks = 1:length(licks)
                        lLickOn(iLicks) = find(botTimes > licks(iLicks),1); %lick onset time
                    end
                    lLickOn = lLickOn(2:end);
                    lLickOn = [lLickOn lLickOn+1];
                end
                
                rLickOn = [];
                if isfield(SessionData.RawEvents.Trial{iFiles}.Events,'Port3In') %check for right licks
                    licks = SessionData.RawEvents.Trial{iFiles}.Events.Port3In;
                    licks(licks < SessionData.RawEvents.Trial{iFiles}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
                    
                    for iLicks = 1:length(licks)
                        rLickOn(iLicks) = find(botTimes > licks(iLicks),1); %lick onset time
                    end
                    rLickOn = rLickOn(2:end);
                    rLickOn = [rLickOn rLickOn+1];
                end
                
                %% collect stim frames and lick frames
                lLickOn(lLickOn > size(mouthTrace,3)) = [];
                rLickOn(rLickOn > size(mouthTrace,3)) = [];
                
                lickFramesL{iFiles} = mouthTrace(:,:,lLickOn);
                lickFramesR{iFiles} = mouthTrace(:,:,rLickOn);
                
                % give some feedback over progress
                if rem(iFiles,50) == 0
                    fprintf(1, 'Current trial is %d out of %d\n', iFiles,length(movieFiles));
                    toc
                end
            end
            
            %% merge video data
            lickFramesL = mat2gray(cat(3,lickFramesL{:}));
            lickFramesR = mat2gray(cat(3,lickFramesR{:}));
            stimFrames = mat2gray(cat(3,stimFrames{:}));
            diffStimFrames = mat2gray(cat(3,diffStimFrames{:}));
            baseFrames = mat2gray(cat(3,baseFrames{:}));
            
            %% create masks and find potential lick events
            fSize = (size(mouthTrace,1) - 1)/2;
            minSize = 200; %minimum tongue size
            erodeSize = 4;
            numSVD = 200;
            
            leftMask = false(size(lickFramesL,1),size(lickFramesL,2));
            leftMask(fSize/5:end-(fSize/5),fSize:(fSize*2)-(fSize/2.5)) = true;

            rightMask = false(size(lickFramesL,1),size(lickFramesL,2));
            rightMask(fSize/5:end-(fSize/5),fSize/2.5:fSize) = true;
            
            fullMask = leftMask | rightMask;
            
            indOut = lickCheck(lickFramesL,leftMask,rightMask,fullMask,minSize,erodeSize);
            leftIdx = indOut > 0;
            
            indOut = lickCheck(lickFramesR,leftMask,rightMask,fullMask,minSize,erodeSize);
            rightIdx = indOut > 0;
            
            indOut= lickCheck(diffStimFrames,leftMask,rightMask,fullMask,minSize,erodeSize);
            stimIdx = indOut > 0;
            
            % make new masks based on average results from left and right licks
            temp = mat2gray(mean(lickFramesL(:,:,leftIdx),3))>0.5;
            a1 = find(sum(temp,1));
            b1 = find(sum(temp,2));

            temp = mat2gray(mean(lickFramesR(:,:,rightIdx),3))>0.5;
            a2 = find(sum(temp,1));
            b2 = find(sum(temp,2));
                        
%             mWidth = min([min([a1 a2]) size(fullMask,2) - max([a1 a2])]); %width of new mask           
%             mWidth = mWidth + 1 : size(fullMask,2) - mWidth;
%             mHeight = min([b1;b2]) : max([b1;b2]);
            mWidth = (size(fullMask,1)-1)/10:size(fullMask,2) - 1-(size(fullMask,1)-1)/10; %width of new mask           
            mHeight = (size(fullMask,1)-1)/2.5 : size(fullMask,1)-1-(size(fullMask,1)-1)/10; %height of new mask
          
            % mask for most of the licking activity. Use this for svd compression.
            fullMask = false(size(fullMask));
            fullMask(mHeight,mWidth) = true;
            
            %collect data and perform svd compression
            lickFramesR = lickFramesR(:,:,rightIdx);
            lickFramesR = cat(3,lickFramesR,flip(lickFramesR,2));
            
            lickFramesL = lickFramesL(:,:,leftIdx);
            lickFramesL = cat(3,lickFramesL,flip(lickFramesL,2));
            
            %%
            data = diffStimFrames(:,:,stimIdx);
            data = cat(3,data,flip(data,2));
            data = cat(3,data,lickFramesL,lickFramesR); %lick data
            data = arrayShrink(data,~fullMask);

            ind = randperm(size(baseFrames,3));
            ind = ind(1:size(data,2));
            data = cat(2,data,arrayShrink(baseFrames(:,:,ind),~fullMask)); %add non-lick examples
            data = cat(2,data,arrayShrink(diffStimFrames,~fullMask)); %add stim period data
            
            % compute U and V
            COV = double(data *  data'/size(data,2));
            [U,Sv] = eigs(COV, opts.nSVD, 'la');
            V = U' * data;
            
            % compute weights
%             [betas, dev, stats] = mnrfit(V(1:100,:)',[ones(length(ind),1);zeros(length(ind),1)]+1);
                        
            svmModel = fitcsvm(V(1:100,1:length(ind)*2)',[true(length(ind),1);false(length(ind),1)],'Standardize',true,'KernelFunction','linear');
            csvmModel = crossval(svmModel);            
            classLoss = kfoldLoss(csvmModel);
            
            
            ind1 = predict(svmModel,V(1:100,(length(ind)*2)+1:end)');
            
            
            
            %%
            
            
            
            
            
            
            lickFramesL = arrayShrink(lickFramesL(:,:,leftIdx),~fullMask);
            lickFramesR = arrayShrink(lickFramesR(:,:,rightIdx),~fullMask);
            
            
            
            
            
            
            leftMask = false(size(lickFramesL,1),size(lickFramesL,2));
            leftMask(min(a) : max(a),min(b) : max(b)) = true;
            
            rightMask = false(size(lickFramesL,1),size(lickFramesL,2));
            rightMask(min(a) : max(a),min(b) : max(b)) = true;
            
            % combine left and right masks and check stim data
            leftMask = leftMask | flip(rightMask,2);
            rightMask = rightMask | flip(leftMask,2);
                       

            
            leftMask(1:end,fSize:(fSize*2)-(fSize/2.5)) = true;
            rightMask(1:end,fSize/2.5:fSize) = true;
            
            
            
            
%             covRatio = (MaskMean+0.0001) ./ (CenterMean+0.0001);


            
            indOut = lickCheck(baseFrames,leftMask,rightMask,fullMask,minSize,erodeSize);
            baseIdx = indOut > 0;            
            
            %% run logitisc regression model

            lWeights = mnrfit(excessRate(trialSelect,ceil(0.05/binSize)+1:end-ceil(0.05/binSize)),lInd(trialSelect)'+1);

            
          %% 
          % change code to only give back regular licks for stim period.
          % left and right licks should come from lick difference movies,
          % so the same analysis can be applied.
          
            %%
            
            if sum(lIdx) > 1000;
                lIdx = indOut == 1 & lRatio < prctile(lRatio(indOut == 1), round((1000 / sum(lIdx)) * 100));
            end
            
            leftData = lickFramesL(:,:,lIdx);
            leftData = reshape(leftData,[],size(leftData,3));
            
            [indOut, rMaskMean, rCenterMean] = lickCheck(lickFramesR,leftMask,rightMask,centerMask,minSize,erodeSize);
            rRatio = (rMaskMean+0.0001) ./ (rCenterMean+0.0001);
            rIdx = indOut == 2 & rRatio < median(rRatio);
            
            rightData = lickFramesR(:,:,rIdx);
            rightData = reshape(rightData,[],size(rightData,3));
            
            [indOut, sMaskMean, sCenterMean] = lickCheck(diffStimFrames,leftMask,rightMask,centerMask,minSize,erodeSize);
            sRatio = (sMaskMean+0.0001) ./ (sCenterMean+0.0001);
            sIdx = indOut == 2 & sRatio < median(sRatio);
            
            %% run classifier based on selected training set

                        stimFrames = reshape(stimFrames,[],size(stimFrames,3));
            
            leftTest = bsxfun(@minus,stimFrames,mean(leftData,2));
            leftTest = bsxfun(@rdivide,leftTest,std(leftData,[],2));
            leftTest = sqrt(sum(leftTest.^2));
            
            rightTest = bsxfun(@minus,stimFrames,mean(rightData));
            rightTest = bsxfun(@rdivide,rightTest,std(rightData));
            rightTest = sqrt(sum(rightTest.^2));
            
            %%
            
            
            
            
            
            %             [leftInd rightInd]= checkLickFrames(stimFrames,leftMask,rightMask);
            
            
            
            %
            %             for iFrames = 1:size(stimFrames,3)
            %                 [~,center(iFrames,:)] = Behavior_EyeCheck(1-stimFrames(:,:,iFrames),0.9,0.5,0); %find areas of interest and get center point
            %             end
            %
            %             % identify frames that have a areas, centered in roughly the right part of the image
            %             leftInd = center(:,1) < fSize*2-(fSize/5) & center(:,1) > fSize & center(:,2) > fSize-(fSize/2.5) & center(:,2) < fSize+(fSize/2.5);
            %             rightInd = center(:,1) > fSize/5 & center(:,1) < fSize & center(:,2) > fSize-(fSize/2.5) & center(:,2) < fSize+(fSize/2.5);
            %
            %             %
            
            
            
            % %% find mask and determine licking threshold
            %             leftMask = false(size(lickFramesL,1),size(lickFramesL,2));
            %             leftMask(fSize-(fSize/2.5):fSize+(fSize/2.5),fSize/5:fSize) = true;
            %             leftThresh = median(mean(arrayShrink(lickFramesL,~leftMask)));
            %             leftIdx = (mean(arrayShrink(stimFrames,~leftMask)) > leftThresh) & (mean(arrayShrink(stimFrames,~leftMask)) > mean(arrayShrink(stimFrames,leftMask)));
            %
            %
            %             rightMask = false(size(lickFramesL,1),size(lickFramesL,2));
            %             rightMask(fSize-(fSize/2.5):fSize+(fSize/2.5),fSize:(fSize*2)-(fSize/5)) = true;
            %             rightThresh = median(mean(arrayShrink(lickFramesL,~rightMask)));
            %             rightIdx = (mean(arrayShrink(stimFrames,~rightMask)) > rightThresh) & (mean(arrayShrink(stimFrames,~rightMask)) > mean(arrayShrink(stimFrames,rightMask)));
            %
            
            %             rightMask = mat2gray(mean(lickFramesR,3));
            %             level = graythresh(rightMask);
            %             rightMask = rightMask < level;
            %             rightMask = bwareaopen(~rightMask,50);
            %             rightMask = imclose(rightMask,strel('disk',4));
            %             rightMask = imerode(rightMask,strel('disk',4));
            %             rightMask = imdilate(rightMask,strel('disk',4));
            %             rightMask = imfill(bwareaopen(rightMask,100),'holes');
            %             rightThresh = median(mean(arrayShrink(lickFramesR,rightMask))); %threshold to count a left lick event
            %
            %             leftMask = mat2gray(mean(lickFramesL,3));
            %             level = graythresh(leftMask);
            %             leftMask = leftMask < level;
            %             leftMask = bwareaopen(~leftMask,50);
            %             leftMask = imclose(leftMask,strel('disk',4));
            %             leftMask = imerode(leftMask,strel('disk',4));
            %             leftMask = imdilate(leftMask,strel('disk',4));
            %             leftThresh = median(mean(arrayShrink(lickFramesL,leftMask))); %threshold to count a left lick event
            
            %%
            
            
            
            %             temp = arrayShrink(lickFramesL,leftMask);
            %
            %
            % cFrame = imadjust(leftMask,[0.5 1],[]); %improve contrast
            %
            %
            %
            %             cFrame = imadjust(cFrame,[0 imThresh],[]); %improve contrast
            % cFrame = cFrame < level;
            % cFrame = imclose(cFrame,strel('disk',4));
            % cFrame = imerode(cFrame,strel('disk',2));
            % cFrame = imdilate(cFrame,strel('disk',2));
            %
            %             a
        end
    end
end

function [indOut, maskOut, centerOut] = lickCheck(framesIn,leftMask,rightMask,centerMask,minSize,erodeSize)
%function to identiy 'lick-like' frames. Frame needs to have a spot of a
%given size in the right place but also some change in the center of the
%image where the tongue originates.

indOut = zeros(1,size(framesIn,3));
maskOut = zeros(1,size(framesIn,3));
centerOut = zeros(1,size(framesIn,3));
level = prctile(framesIn(:),90);

for iFrames = 1:size(framesIn,3)
% for iFrames = 1:1000
    
    %     cFrame = imadjust(framesIn(:,:,iFrames),[0 0.5],[]); %adjust contrast
    cFrame = framesIn(:,:,iFrames);
%     subplot(1,2,1);
%     imagesc(framesIn(:,:,iFrames)); colormap jet; axis image;
    
    cFrame =cFrame > prctile(cFrame(:),90);
    cFrame = imerode(cFrame,strel('disk',erodeSize));
    cFrame = bwareaopen(imfill(cFrame,'holes'),minSize);
    
    areaInfo = regionprops( cFrame, 'Centroid','MajorAxisLength','MinorAxisLength','Solidity');
    
    if length(areaInfo) == 1
%         cntr = .5 * [size(cFrame,2) size(cFrame,1)]; % X-Y coordinates and NOT Row/Col
%         d = sqrt( sum( bsxfun(@minus,vertcat(areaInfo.Centroid), cntr ).^2, 2 ) );
%         [~, idx] = min(d);
        idx = 1;
        
        center = round(areaInfo(idx).Centroid);
        axes = [areaInfo(idx).MinorAxisLength areaInfo(idx).MajorAxisLength];
        if leftMask(center(2),center(1)) && areaInfo(idx).Solidity > .75 && min(axes/max(axes)) > 0.5
            indOut(iFrames) = 1;
            maskOut(iFrames) = mean(arrayShrink(framesIn(:,:,iFrames),~leftMask));
            centerOut(iFrames) = mean(arrayShrink(framesIn(:,:,iFrames),leftMask));
            
        elseif rightMask(center(2),center(1)) && areaInfo(idx).Solidity > .75 && min(axes/max(axes)) > 0.5
            indOut(iFrames) = 2;
            maskOut(iFrames) = mean(arrayShrink(framesIn(:,:,iFrames),~rightMask));
            centerOut(iFrames) = mean(arrayShrink(framesIn(:,:,iFrames),rightMask));
        end
    end
    
%     subplot(1,2,2);
%     imagesc(cFrame); colormap jet; axis image;
%     disp([num2str(iFrames) ' - ' num2str(indOut(iFrames)) ]);
%     pause
end
