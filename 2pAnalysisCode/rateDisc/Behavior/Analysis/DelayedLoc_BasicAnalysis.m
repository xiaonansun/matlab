function [Performance,bhv] = DelayedLoc_BasicAnalysis(Animal,cPath,lSessions,binSize,highDetection,showPlot,singleMod,minDelay)
% Analze behavioral data from delayed localization paradigm to get basic readout of animal performance.
% Reads all available files in 'cPath', ignoring data files that contain less then 100 trials.
%
% Optional inputs:
% lSessions: Only use last 'lSessions' sessions for behavioral analysis.
% binSize: Number of bins used to compute delay performance curves. The default is 10.
% highDetection: Only use sessions at which detection was at 90% or higher.
% singleMod: Only use sessions at which stimuli of stimType == 'singleMod' were presented. The default is all trials.
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

if ~exist('singleMod','var') || isempty(singleMod)
    singleMod = 0;
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
decWait = 1; %maximum decision time, used to compute session performance

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
earlyStimTime = 0.5; %separator to count events in early part of the stimulus. events after earlyStimTime are counted for the late part.

%% load data
bhv = []; Cnt = 0;

for iFiles = 1:size(Files,1)
    
    load([cPath Files(iFiles,:)], 'SessionData'); %load current bhv file
    if isfield(SessionData,'Rewarded')
        SessionData.Rewarded = logical(SessionData.Rewarded);
    end
    
    useData = isfield(SessionData,'decisionGap') && length(SessionData.Rewarded) > 100; % if file contains at least 100 trials
    if singleMod > 0
        useData = isfield(SessionData,'decisionGap') && length(SessionData.Rewarded) > 100 && ~any(SessionData.StimType ~= singleMod); % if file contains at least 100 trials and trials of type singleMod
    end
    
    if minDelay > 0
        useData = isfield(SessionData,'decisionGap') && length(SessionData.Rewarded) > 100 && any(SessionData.decisionGap > minDelay); % if file contains at least 100 trials and sufficiently long waiting periods
    end
    
    if useData
        Cnt = Cnt+1;
        DayNr(Cnt) = cDate(iFiles);

        %% get some single session performance data
        ind = logical(SessionData.Assisted) & SessionData.decisionGap <= decWait; %index for trials that were self-performed
        lInd = SessionData.CorrectSide == 1; %index for left-choice trials
        Performance.SelfPerformed(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        Performance.AllTrials(1,Cnt) = sum(SessionData.Rewarded)/sum(SessionData.Rewarded+SessionData.Punished); %performance in all trials
        Performance.LeftPerformed(1,Cnt) = sum(SessionData.Rewarded(ind & lInd))/sum(SessionData.Rewarded(ind & lInd)+SessionData.Punished(ind & lInd));
        Performance.RightPerformed(1,Cnt) = sum(SessionData.Rewarded(ind & ~lInd))/sum(SessionData.Rewarded(ind & ~lInd)+SessionData.Punished(ind & ~lInd));
        Performance.Date{1,Cnt} = datestr(cDate(Cnt));
        
        ind =  logical(SessionData.Assisted) & SessionData.decisionGap <= decWait; %index for detection trials
        Performance.Detection(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        Performance.LeftDetected(1,Cnt) = sum(SessionData.Rewarded(ind & lInd))/sum(SessionData.Rewarded(ind & lInd)+SessionData.Punished(ind & lInd));
        Performance.RightDetected(1,Cnt) = sum(SessionData.Rewarded(ind & ~lInd))/sum(SessionData.Rewarded(ind & ~lInd)+SessionData.Punished(ind & ~lInd));
        
        ind =  logical(SessionData.Assisted) & SessionData.StimType == 1 & SessionData.decisionGap <= decWait; %index for vision trials
        Performance.Vision(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        
        ind =  logical(SessionData.Assisted) & SessionData.StimType == 2 & SessionData.decisionGap <= decWait; %index for audio trials
        Performance.Audio(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        
        ind =  logical(SessionData.Assisted) & SessionData.StimType == 3 & SessionData.decisionGap <= decWait; %index for multisensory trials
        Performance.Mixed(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials        

        %% compute reaction times
%         for iTrials = 1:SessionData.nTrials
%             if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
%                 lLick = SessionData.RawEvents.Trial{iTrials}.Events.Port1In; %left licks
%             else
%                 lLick = NaN; %no licks
%             end
%             if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port3In') %check for licks
%                 rLick = SessionData.RawEvents.Trial{iTrials}.Events.Port3In; %right licks
%             else
%                 rLick = NaN; %no licks
%             end
%             
%             if (any(isnan(lLick)) && any(isnan(rLick))) || any(any(isnan(SessionData.RawEvents.Trial{iTrials}.States.WaitForResponse)))
%                 reactionTime(iTrials) = NaN; %no response
%             else
%                 window = SessionData.RawEvents.Trial{iTrials}.States.WaitForResponse;
%                 lLick = lLick-window(1); lLick(lLick < 0 | lLick > (window(2) - window(1))) = [];lLick(isnan(lLick)) = inf; %compute left response during decision window
%                 rLick = rLick-window(1); rLick(rLick < 0 | rLick > (window(2) - window(1))) = [];rLick(isnan(rLick)) = inf; %compute right response during decision window
%                 
%                 if min(lLick) <  min(rLick)  % if left licks are first
%                     reactionTime(iTrials) = lLick(1);
%                 elseif  min(lLick) >  min(rLick)   % if right licks are first
%                     reactionTime(iTrials) = rLick(1);
%                 else
%                     reactionTime(iTrials) = NaN;
%                 end
%                 clear window
%             end
%             clear  lLick rLick
%         end
%         
%         ind = logical(SessionData.Assisted); %index for trials that were self-performed
%         Performance.rTime(1,Cnt).All = nanmean(reactionTime); %mean reaction time for session
%         Performance.rTime(1,Cnt).Correct = nanmean(reactionTime(SessionData.Rewarded(ind))); %mean reaction time in correct trials
%         Performance.rTime(1,Cnt).False = nanmean(reactionTime(SessionData.Punished(ind))); %mean reaction time in correct trials
        
        %% combine into one larger array
%         SessionData.rTime = reactionTime; clear reactionTime
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

% if length(sessionSelect) < lSessions
%     lSessions = length(sessionSelect);
% end
% sessionSelect = sessionSelect(1:lSessions); %only use last 'lSessions'
% disp(['Selected ' num2str(lSessions) '/' num2str(Cnt) ' remaining sessions.']);

%only use last 'lSessions' days
bhv = selectBehavior(bhv,sessionSelect); %only use trials from selecteded sessions
Performance = selectBehavior(Performance,sessionSelect); %only use performance from selecteded sessions

DayNr = DayNr(sessionSelect) - DayNr(end) + 1; %convert Days into relative values, starting at 1 for first training day
SessionNr = length(DayNr):-1:1;
fDay = Performance.Date{end}; %first date in dataset
lDay = Performance.Date{1}; %last date in dataset
disp(['First date: ' fDay]);
disp(['Last date: ' lDay]);

%% compute performance for different decision delays
bhv = selectBehaviorTrials(bhv,bhv.decisionGap <= maxDelay); %only use trials that are below maximum allowed delay

if length(unique(bhv.decisionGap)) <= 10
    gaps = unique(bhv.decisionGap(bhv.decisionGap > 0));
else
    [~,gaps] = histcounts((bhv.decisionGap(bhv.decisionGap > 0)),binSize);
    gaps(end) = [];
end
gapStep = mean(diff(gaps));
if isnan(gapStep) || isinf(gapStep)
    gapStep = 1;
end

% compute zero delay detection
ind = bhv.decisionGap == 0 & bhv.Assisted; %find trials with right gap values
ind1 = ind & bhv.StimType == 1; %vision
ind2 = ind & bhv.StimType == 2; %audio
ind3 = ind & bhv.StimType == 3; %audiovisual
allInd = ind1 | ind2 | ind3;

% no delay
[discPerf{1}(1,1),distConv{1}(:,1,1),Performance.tCount{1}(1,1)] = computeBehavior(bhv,maxTrialCnt,allInd); %compute performance and error for all trials
[discPerf{1}(1,2),distConv{1}(:,1,2),Performance.tCount{1}(1,2)] = computeBehavior(bhv,maxTrialCnt,allInd & bhv.CorrectSide == 1); %compute performance and error for left trials
[discPerf{1}(1,3),distConv{1}(:,1,3),Performance.tCount{1}(1,3)] = computeBehavior(bhv,maxTrialCnt,allInd & bhv.CorrectSide == 2); %compute performance and error for right trials

% isolate different modalites
[Performance.modDiscPerf{1}(1,1),modDistConv{1}(:,1,1),Performance.modtCount{1}(1,1)] = computeBehavior(bhv,maxTrialCnt,ind1); %vision trials
[Performance.modDiscPerf{1}(1,2),modDistConv{1}(:,1,2),Performance.modtCount{1}(1,2)] = computeBehavior(bhv,maxTrialCnt,ind2); %audio trials
[Performance.modDiscPerf{1}(1,3),modDistConv{1}(:,1,3),Performance.modtCount{1}(1,3)] = computeBehavior(bhv,maxTrialCnt,ind3); %mixed trials

for iGaps = 1:length(gaps)
    
    ind = bhv.decisionGap >= gaps(iGaps) & bhv.decisionGap < (gaps(iGaps) +gapStep) & bhv.Assisted; %find trials with right gap values
    ind1 = ind & bhv.StimType == 1; %vision
    ind2 = ind & bhv.StimType == 2; %audio
    ind3 = ind & bhv.StimType == 3; %audiovisual
    allInd = ind1 | ind2 | ind3;
    
    if sum(ind1) < 10; ind1 = false(1,length(ind)); end
    if sum(ind2) < 10; ind2 = false(1,length(ind)); end
    if sum(ind3) < 10; ind3 = false(1,length(ind)); end
    
    [discPerf{1}(iGaps+1,1),distConv{1}(:,iGaps+1,1),Performance.tCount{1}(iGaps+1,1)] = computeBehavior(bhv,maxTrialCnt,allInd); %compute performance and error for all trials
    [discPerf{1}(iGaps+1,2),distConv{1}(:,iGaps+1,2),Performance.tCount{1}(iGaps+1,2)] = computeBehavior(bhv,maxTrialCnt,allInd & bhv.CorrectSide == 1); %compute performance and error for left trials
    [discPerf{1}(iGaps+1,3),distConv{1}(:,iGaps+1,3),Performance.tCount{1}(iGaps+1,3)] = computeBehavior(bhv,maxTrialCnt,allInd & bhv.CorrectSide == 2); %compute performance and error for right trials
    
    % isolate different modalites
    [Performance.modDiscPerf{1}(iGaps+1,1),modDistConv{1}(:,iGaps+1,1),Performance.modtCount{1}(iGaps+1,1)] = computeBehavior(bhv,maxTrialCnt,ind1); %vision trials
    [Performance.modDiscPerf{1}(iGaps+1,2),modDistConv{1}(:,iGaps+1,2),Performance.modtCount{1}(iGaps+1,2)] = computeBehavior(bhv,maxTrialCnt,ind2); %audio trials
    [Performance.modDiscPerf{1}(iGaps+1,3),modDistConv{1}(:,iGaps+1,3),Performance.modtCount{1}(iGaps+1,3)] = computeBehavior(bhv,maxTrialCnt,ind3); %mixed trials

end

if length(unique(bhv.decisionGap)) > 10
    gaps = gaps +gapStep/2;
end
gaps = [0 gaps];
Performance.gaps = gaps;

if showPlot
%% Overview figure for sessions
figure('name',[Animal ' - Learning curves; Start date: ' fDay ' ; End date: ' lDay])
subplot(1,2,1)
Data = Performance.Detection;
plot(SessionNr,Data,'-ok','linewidth',2);hold on
plot(SessionNr(~isnan(Data)),smooth(Data(~isnan(Data))),'--k','linewidth',2)

line([0 length(SessionNr)+1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
title([Animal ' - Detection Performance']); xlabel('#sessions','fontsize',15); ylabel('performance (%)','fontsize',15);
axis square; set(gca,'FontSize',12); ylim([0.4 1.05]);
% legend({'Combined'},'Location','SouthEast')

subplot(1,2,2)
Data = Performance.Mixed;
plot(SessionNr,Data,'-ok','linewidth',2);hold on
Data = Performance.Vision;
plot(SessionNr,Data,'-og','linewidth',2);hold on
Data = Performance.Audio;
plot(SessionNr,Data,'-or','linewidth',2);hold on

line([0 length(SessionNr)+1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
title([Animal ' - Detection Performance']); xlabel('#sessions','fontsize',15); ylabel('performance (%)','fontsize',15);
axis square; ylim([0.35 1]); set(gca,'FontSize',12);
legend({'Multisensory' 'Vision' 'Audio'})
% legend({'All' 'Left' 'Right'},'Location','SouthEast')

%% overview for distractor performance
h = figure('name',[Animal ' - Delay curves; Start date: ' fDay ' ; End date: ' lDay]);
ax = cla(h);hold on;
axis square; ylim([0.35 1]); xlim([gaps(1)-gapStep gaps(end)+gapStep]); set(gca,'FontSize',12);
line(xlim(ax),repmat(discPerf{1}(1),1,2),'linestyle','-','linewidth',2,'color',[0.5 0.5 0.5])
line(xlim(ax),[0.5 0.5],'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5])
title([Animal ' - Discrimination Performance']); xlabel('decision delay (s)','fontsize',15);
ylabel('performance (%)','fontsize',15);
    
plotDist = gaps(~isnan(discPerf{1}(:,1)));
plotPerf = discPerf{1}(~isnan(discPerf{1}(:,1)),:);
plotError = distConv{1}(:,~isnan(discPerf{1}(:,1)),:);

if ~isempty(plotPerf)
    if size(plotError,2) == 1
        plotError(1,:,:) = plotPerf - squeeze(plotError(1,:,:))';
        plotError(2,:,:) = squeeze(plotError(2,:,:))' - plotPerf;
    else
        plotError(1,:,:) = plotPerf - squeeze(plotError(1,:,:));
        plotError(2,:,:) = squeeze(plotError(2,:,:)) - plotPerf;
    end
    
    errorbar(plotDist,plotPerf(:,2),plotError(1,:,2),plotError(2,:,2), ...
        '--og','Markersize',5,'MarkerEdgeColor','g','MarkerFaceColor','w','linewidth',2)
    
    errorbar(plotDist,plotPerf(:,3),plotError(1,:,3),plotError(2,:,3), ...
        '--or','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','w','linewidth',2)
    
    errorbar(plotDist,plotPerf(:,1),plotError(1,:,1),plotError(2,:,1), ...
        '-ok','Markersize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2);

end


%% compare performance at different modalities
if ~isempty(Performance.modDiscPerf{1}(~isnan(Performance.modDiscPerf{1}(:,1:2)))) %single modality trials in dataset
    
    h = figure;
    title([Animal ' - Single modalities']);
    ax = cla(h);hold on;
    axis square; ylim([0.35 1]); xlim([gaps(1)-gapStep gaps(end)+gapStep]); set(gca,'FontSize',12);
    line(xlim(ax),[0.5 0.5],'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5])
    title([Animal ' - Discrimination Performance']); xlabel('decision delay (s)','fontsize',15);
    ylabel('performance (%)','fontsize',15);
    
    plotError = modDistConv{1};
    plotDist = gaps;
    plotPerf = Performance.modDiscPerf{1};
    
    if size(plotError,2) == 1
        plotError(1,:,:) = plotPerf - squeeze(plotError(1,:,:))';
        plotError(2,:,:) = squeeze(plotError(2,:,:))' - plotPerf;
    else
        plotError(1,:,:) = plotPerf - squeeze(plotError(1,:,:));
        plotError(2,:,:) = squeeze(plotError(2,:,:)) - plotPerf;
    end
    
    errorbar(plotDist,plotPerf(:,1),plotError(1,:,1),plotError(2,:,1), ...
        '-og','Markersize',9,'MarkerEdgeColor','g','MarkerFaceColor','w','linewidth',2)
    
    errorbar(plotDist,plotPerf(:,2),plotError(1,:,2),plotError(2,:,2), ...
        '-or','Markersize',9,'MarkerEdgeColor','r','MarkerFaceColor','w','linewidth',2)
    
    errorbar(plotDist,plotPerf(:,3),plotError(1,:,3),plotError(2,:,3), ...
        '-ok','Markersize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2);  
    
end
end
end