function Behavior_CheckVideoTiming(Animal,Rec)

Paradigm = 'SpatialDisc';
cFolder = 'BehaviorVideo';
cPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep cFolder filesep]; %Widefield data path

%% load bpod file
bhvFile = dir([fileparts(cPath(1:end-1)) filesep Animal '*' Paradigm '*.mat']);
load([fileparts(cPath(1:end-1)) filesep bhvFile.name])

%% find nr of cams
files = dir([cPath '*frameTimes*.mat']);
temp = cat(1,files.name);
ind = strfind(temp(1,:),'_');ind = ind(end);
nrCams = length(unique(str2num(temp(:,ind+1))));

%%check for LED indicator and ask for location if required
for iCams = 1:nrCams
if ~isfield(SessionData,['led' int2str(iCams) 'Pos'])
    movieFiles = dir([cPath '*Video*_' int2str(iCams) '.mj2']); %all files for current cam based on integer before .mj2
    
    checker = true;
    while checker
        rawData = squeeze(importdata([cPath movieFiles(randperm(length(movieFiles),1)).name]));
        [~,~,~,~,tracePosition] = showStack(rawData);
        Wait = input('Happy with indicator position? - Type "y" to continue\n','S');
        if strcmpi(Wait,'y')
            checker = false;
            tracePosition(1:2) = tracePosition(1:2) + tracePosition(3)/2;
            SessionData.(['led' int2str(iCams) 'Pos']) = tracePosition(1:3);
        end
    end
end
end
save([fileparts(cPath(1:end-1)) filesep bhvFile.name],'SessionData'); %save bhv file with indicator position

%% loop through available webcams
for iCams = 1:nrCams
    
    timeFiles = dir([cPath '*frameTimes*_' int2str(iCams) '.mat']); %all files for current cam based on integer before .mat
    movieFiles = dir([cPath '*Video*_' int2str(iCams) '.mj2']); %all files for current cam based on integer before .mj2
    ledPos = SessionData.(['led' int2str(iCams) 'Pos']);
    Cnt = 0;
    
    for iFiles = 1:length(movieFiles)
        tic
        load([cPath timeFiles(iFiles).name]); %load frameTimes
        frameTimes = (frameTimes - SessionData.TrialStartTime(iFiles))  * (86400); %convert to s, relative to trialstart

        rawData = squeeze(importdata([cPath movieFiles(iFiles).name])); %get video data and isolate indicator trace
        ledTrace = rawData((ledPos(2)-ledPos(3)-1) + (1:ledPos(3)*2), (ledPos(1)-ledPos(3)-1) + (1:ledPos(3)*2), :);
        ledTrace = mean(reshape(ledTrace,(ledPos(3)*2)^2,[]));
        indOn = zscore(ledTrace) > 1;
        
        onsets = SessionData.RawEvents.Trial{iFiles}.States.Reset(:,1);
        for iInds = 1:length(onsets)
            temp = frameTimes - onsets(iInds);
            trigOn{iCams}(Cnt + iInds) = find(temp > 0,1); %expected frame where trigger should be first visible
            trigTimeDiff{iCams}(Cnt + iInds) = temp(trigOn{iCams}(Cnt + iInds)); %time difference between indicator and next acquired frame
           
            temp = find(indOn(trigOn{iCams}(Cnt + iInds) : end), 1) - 1;
            if isempty(temp)
                trigTimeShift{iCams}(Cnt + iInds) =NaN;
            else
                trigTimeShift{iCams}(Cnt + iInds) = temp; %frame shift between indicator and observed event
            end
        end
        Cnt = Cnt + iInds;
        toc
    end
end
save([fileparts(cPath(1:end-1)) filesep 'camTiming'],'trigOn','trigTimeDiff','trigTimeShift');

%% make figures
figure
for x = 1:2
subplot(1,2,x);
histogram(trigTimeDiff{x}(trigTimeShift{x} == 0),10); hold
histogram(trigTimeDiff{x}(trigTimeShift{x} == 1),10)
histogram(trigTimeDiff{x}(trigTimeShift{x} == 2),10)
legend({'0' '1' '2'});
axis square; title(['Cam ' int2str(x) ' - high delay: ' int2str(sum((trigTimeShift{1} > 2 | isnan(trigTimeShift{1})))) '\' int2str(length(trigTimeShift{1}))]);
xlim([-0.0025 0.0425])
end