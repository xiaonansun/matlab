function [Performance,bhv] = DelayedLoc_learningCurves(Animal,cPath,lSessions,binSize,highDetection,showPlot,minDelay)
% Analze behavioral data from delayed localization paradigm to get basic readout of animals discrimination performance.
% Reads all available files in 'cPath', ignoring data files that contain less then 100 trials.
%
% Optional inputs:
% lSessions: Only use last 'lSessions' sessions for behavioral analysis.
% binSize: Number of bins used to compute discrimination. Default is 10.
% highDetection: Only use sessions at which detection was at 90% or higher.
% minDelay: Only use sessions at which at least one decisionGap > 'minDelay' was presented. The default is all trials.
% showPlot: Plot behavioral results. Default is true.

%% check optional input
if ~exist('lSessions','var') || isempty(lSessions)
    lSessions = inf;
end

if ~exist('highDetection','var') || isempty(highDetection)
    highDetection = 0;
end

if ~exist('showPlot','var') || isempty(showPlot)
    showPlot = true;
end

if ~exist('minDelay','var') || isempty(minDelay)
    minDelay = 0;
end

if ~exist('binSize','var') || isempty(binSize)
    binSize = 10;
end

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end
modId = [2 4 6]; %id for modalities. Default is audio, somatosensory and audiosomatosensory
modLabels = {'audio' 'tactile' 'audiotactile'}; %labels for different modalities
modColor = ['r' 'b' 'k']; %colors for different modalities
nrBlocks = 5; %number of blocks the data is separated in later

%% assign some variables
paradigm = 'SpatialDisc';
for iChecks = 1:10 %check for files repeatedly. Sometimes the server needs a moment to be indexed correctly
    Files = ls([cPath '\' Animal '\' paradigm '\Session Data\*' Animal '*']); %behavioral files in correct cPath
    if ~isempty(Files)
        break;
    end
    pause(0.1);
end

cPath = [cPath Animal '\' paradigm '\Session Data\']; %folder with behavioral data
maxTrialCnt = 1000; %maximum trials per datapoint
maxDelay = 3; %maximum delay in seconds

cDate = datenum(Files(:,length([Animal '_' paradigm '_'])+1:length([Animal '_' paradigm '_'])+10)); %isolate dates from Filenames
cDate = cDate + (str2num(Files(:,length([Animal '_' paradigm '_'])+19:end-4))*0.01); %add session nr to timestamp

[cDate,ind] = sort(cDate,'descend'); %sort counts to get the order of files to days correct. Newest file should be first in the list.
cDate = floor(cDate);
Files = Files(ind,:); %adjust order of filenames to get it to be chronological

%% load data
bhv = []; Cnt = 0;

for iFiles = 1:size(Files,1)
    try
    load([cPath Files(iFiles,:)], 'SessionData'); %load current bhv file
    catch
    end
    
    if isfield(SessionData,'Rewarded')
        SessionData.Rewarded = logical(SessionData.Rewarded);
    end
    
    useData = isfield(SessionData,'Rewarded') && length(SessionData.Rewarded) > 100; % if file contains at least 100 trials and different distractor rates
    
    if minDelay > 0
        useData = isfield(SessionData,'decisionGap') && length(SessionData.Rewarded) > 100 && any(SessionData.decisionGap > minDelay); % if file contains at least 100 trials and sufficiently long waiting periods
    end
    
    
    if useData
        Cnt = Cnt+1;
        DayNr(Cnt) = cDate(iFiles);
        Performance.Date{1,Cnt} = datestr(cDate(iFiles));
            
        for iMod = 1:3
            %% get some single session performance data
            ind = ~SessionData.DidNotChoose & taertaert~SessionData.DidNotLever & logical(SessionData.Assisted) & SessionData.StimType == modId(iMod); %only use active trials
            Performance.SelfPerformed(iMod,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
            
            dInd = logical(SessionData.Assisted & SessionData.DistStim == 0); %index for trials that were detection only
            Performance.Detection(iMod,Cnt) = sum(SessionData.Rewarded(ind & dInd))/sum(SessionData.Rewarded(ind & dInd)+SessionData.Punished(ind & dInd)); %peformance for detection trials
            Performance.Discrimination(iMod,Cnt) = sum(SessionData.Rewarded(ind & ~dInd))/sum(SessionData.Rewarded(ind & ~dInd)+SessionData.Punished(ind & ~dInd)); %peformance for discrimination trials
            Performance.DetTrials(iMod,Cnt) = sum(ind & dInd); %number of detection trials
            Performance.DiscTrials(iMod,Cnt) = sum(ind & ~dInd); %number of discrimination trials
            Performance.AllTrials(Cnt) = length(SessionData.Rewarded);
            
            lInd = SessionData.CorrectSide == 1; %index for left-choice trials
            Performance.LeftPerformed(iMod,Cnt) = sum(SessionData.Rewarded(ind & lInd))/sum(SessionData.Rewarded(ind & lInd)+SessionData.Punished(ind & lInd));
            Performance.RightPerformed(iMod,Cnt) = sum(SessionData.Rewarded(ind & ~lInd))/sum(SessionData.Rewarded(ind & ~lInd)+SessionData.Punished(ind & ~lInd));
        end
        
        %% combine into one larger array
        SessionData.SessionNr = repmat(Cnt,1,SessionData.nTrials); %tag all trials in current dataset with session nr
        bhv = appendBehavior(bhv,SessionData); %append into larger array
        
        if Cnt >= lSessions
            break;
        end
    end
end
disp(['Current subject: ' Animal '; Using ' num2str(Cnt) '/' num2str(size(Files,1)) ' files']);

if Cnt > 0
    %% check for last and high performance sessions
    sessionSelect = 1:Cnt; %all sessions
    if highDetection > 0
        lowInd = Performance.Detection < highDetection; %find sessions with low detection performance
    else
        lowInd = false(1,Cnt);
    end
    sessionSelect(lowInd) = []; %don't use low performance sessions of selected
    disp(['Rejected ' num2str(sum(lowInd)) '/' num2str(length(lowInd))  ' files for detection performance below ' num2str(highDetection*100) '%.']);
    
    %only use last 'lSessions' days
    bhv = selectBehavior(bhv,sessionSelect); %only use trials from selecteded sessions
    Performance = selectBehavior(Performance,sessionSelect); %only use performance from selecteded sessions
    
    DayNr = DayNr(sessionSelect) - DayNr(end) + 1; %convert Days into relative values, starting at 1 for first training day
    SessionNr = length(DayNr):-1:1;
    Performance.SessionNr = SessionNr;
    fDay = Performance.Date{end}; %first date in dataset
    lDay = Performance.Date{1}; %last date in dataset
    disp(['First date: ' fDay]);
    disp(['Last date: ' lDay]);
    
    %% make plots
    if showPlot
        %% cross session - detection / discrimination
        figure('name',[Animal ' - Learning curves; Start date: ' fDay ' ; End date: ' lDay])
        subplot(2,1,1); hold on
        for iMod = 1:3
            plot(SessionNr,Performance.Detection(iMod,:),'-o','linewidth',2, 'color', modColor(iMod));hold on
        end
        
        line([0 length(SessionNr)+1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
        title([Animal ' - Detection / Sessions']); xlabel('#sessions','fontsize',15); ylabel('performance (%)','fontsize',15);
        axis square; set(gca,'FontSize',12); ylim([0.4 1.05]);
        xlim([0 size(Performance.Detection,2)+1])
        
        subplot(2,1,2); hold on
        for iMod = 1:3
            plot(SessionNr,Performance.Discrimination(iMod,:),'--o','linewidth',2, 'color', modColor(iMod));
        end
        line([0 length(SessionNr)+1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
        title([Animal ' - Discrimination / Sessions']); xlabel('#sessions','fontsize',15); ylabel('performance (%)','fontsize',15);
        axis square; set(gca,'FontSize',12); ylim([0.4 1.05]);
        xlim([0 size(Performance.Detection,2)+1])
        
     
    end
end
end