function [Performance,bhv] = SpatialDisc_BasicAnalysis(Animal,path,lSessions,highDetection,singleMod)
% Analze behavioral data from SpatialDisc paradigm to get basic readout of animal performance.
% Reads all available files in 'path', ignoring data files that contain less then 100 trials.
%
% Optional inputs:
% lSessions: Only use last 'lSessions' sessions for behavioral analysis.
% highDetection: Only use sessions at which detection was at 90% or higher.
% singleMod: Only use sessions at which stimuli of stimType == 'singleMod' were presented. The default is all trials.

%% check optional input
if ~exist('lSessions','var') || isempty(lSessions)
    lSessions = inf;
end

if ~exist('highDetection','var') || isempty(highDetection)
    highDetection = 0;
end

if ~exist('singleMod','var') || isempty(singleMod)
    singleMod = 0;
end

%% assign some variables
paradigm = 'SpatialDisc';
Files = ls([path '\' Animal '\' paradigm '\Session Data\*' Animal '*']); %behavioral files in correct path
path = [path Animal '\' paradigm '\Session Data\']; %folder with behavioral data
maxTrialCnt = 1000; %maximum trials per datapoint

cDate = datenum(Files(:,length([Animal '_' paradigm '_'])+1:length([Animal '_' paradigm '_'])+10)); %isolate dates from Filenames
cDate = cDate + (str2num(Files(:,length([Animal '_' paradigm '_'])+19:end-4))*0.01); %add session nr to timestamp

[cDate,ind] = sort(cDate,'descend'); %sort counts to get the order of files to days correct. Newest file should be first in the list.
cDate = floor(cDate);
Files = Files(ind,:); %adjust order of filenames to get it to be chronological
earlyStimTime = 0.5; %separator to count events in early part of the stimulus. events after earlyStimTime are counted for the late part.

%% load data
bhv = []; Cnt = 0;

for iFiles = 1:size(Files,1)
    
    load([path Files(iFiles,:)]); %load current bhv file
    if  isfield(SessionData,'Rewarded')
        SessionData.Rewarded = logical(SessionData.Rewarded);
    end
    
    if isfield(SessionData,'decisionGap') && any(SessionData.decisionGap > 0)
        useData = false; %don't use session if there are any decision gaps
    else
        if singleMod > 0
            useData = isfield(SessionData,'Rewarded') && length(SessionData.Rewarded) > 100 && ~any(SessionData.StimType ~= singleMod); % if file contains at least 100 trials and trials of type singleMod
        else
            useData = isfield(SessionData,'Rewarded') && length(SessionData.Rewarded) > 100; % if file contains at least 100 trials
        end
    end
    
    if useData
        Cnt = Cnt+1;
        DayNr(Cnt) = cDate(iFiles);

        %% get some single session performance data
        ind = logical(SessionData.Assisted); %index for trials that were self-performed
        lInd = SessionData.CorrectSide == 1; %index for left-choice trials
        Performance.SelfPerformed(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        Performance.AllTrials(1,Cnt) = sum(SessionData.Rewarded)/sum(SessionData.Rewarded+SessionData.Punished); %performance in all trials
        Performance.Notes{1,Cnt} = SessionData.Notes{1}; %only use first note if there was anything intersting going on
        Performance.LeftPerformed(1,Cnt) = sum(SessionData.Rewarded(ind & lInd))/sum(SessionData.Rewarded(ind & lInd)+SessionData.Punished(ind & lInd));
        Performance.RightPerformed(1,Cnt) = sum(SessionData.Rewarded(ind & ~lInd))/sum(SessionData.Rewarded(ind & ~lInd)+SessionData.Punished(ind & ~lInd));
        Performance.Date{1,Cnt} = datestr(cDate(Cnt));
        
        ind =  logical(SessionData.Assisted) & SessionData.DistStim == 0; %index for detection trials
        Performance.Detection(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        Performance.LeftDetected(1,Cnt) = sum(SessionData.Rewarded(ind & lInd))/sum(SessionData.Rewarded(ind & lInd)+SessionData.Punished(ind & lInd));
        Performance.RightDetected(1,Cnt) = sum(SessionData.Rewarded(ind & ~lInd))/sum(SessionData.Rewarded(ind & ~lInd)+SessionData.Punished(ind & ~lInd));
        
        ind =  logical(SessionData.Assisted) & SessionData.StimType == 1 & SessionData.DistStim == 0; %index for vision trials
        Performance.Vision(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        
        ind =  logical(SessionData.Assisted) & SessionData.StimType == 2 & SessionData.DistStim == 0; %index for audio trials
        Performance.Audio(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials
        
        ind =  logical(SessionData.Assisted) & SessionData.StimType == 3 & SessionData.DistStim == 0; %index for multisensory trials
        Performance.Mixed(1,Cnt) = sum(SessionData.Rewarded(ind))/sum(SessionData.Rewarded(ind)+SessionData.Punished(ind)); %peformance for self-performed trials        
        
        %% check for occurence of stimulus events in each trial
        SessionData.earlyTarget = cell(1,size(SessionData.Rewarded,2));
        SessionData.lateTarget = cell(1,size(SessionData.Rewarded,2));
        SessionData.earlyDist = cell(1,size(SessionData.Rewarded,2));
        SessionData.lateDist = cell(1,size(SessionData.Rewarded,2));
        
        for iTrials = 1:size(SessionData.Rewarded,2)
            tCnt = 0;
            for x = 1:2:5
                tCnt = tCnt +1;
                if size(SessionData.stimEvents{iTrials},2) >= x
                    
                    [~, b] = max([length(SessionData.stimEvents{iTrials}{x}) length(SessionData.stimEvents{iTrials}{x+1})]);
                    
                    SessionData.allTarget(tCnt,iTrials) = max([length(SessionData.stimEvents{iTrials}{x}) length(SessionData.stimEvents{iTrials}{x+1})]);
                    SessionData.allDist(tCnt,iTrials) = min([length(SessionData.stimEvents{iTrials}{x}) length(SessionData.stimEvents{iTrials}{x+1})]);
                    SessionData.allDisc(tCnt,iTrials) =  SessionData.allDist(tCnt,iTrials) / SessionData.allTarget(tCnt,iTrials);
                    
                    temp1 = sum(SessionData.stimEvents{iTrials}{x} <= earlyStimTime);
                    temp2 = sum(SessionData.stimEvents{iTrials}{x + 1} <= earlyStimTime);
                    
                    temp3 = sum(SessionData.stimEvents{iTrials}{x} > earlyStimTime);
                    temp4 = sum(SessionData.stimEvents{iTrials}{x + 1} > earlyStimTime);
                    
                    if b == 1
                        SessionData.earlyDisc(tCnt,iTrials) = temp2 / temp1;
                        SessionData.lateDisc(tCnt,iTrials) = temp4 / temp3;
                    else
                        SessionData.earlyDisc(tCnt,iTrials) = temp1 / temp2;
                        SessionData.lateDisc(tCnt,iTrials) = temp3 / temp4;
                    end
                
                else
                    SessionData.allTarget(tCnt,iTrials) = NaN;
                    SessionData.allDist(tCnt,iTrials) = NaN;
                    SessionData.allDisc(tCnt,iTrials) = NaN;
                    SessionData.earlyDisc(tCnt,iTrials) = NaN;
                    SessionData.lateDisc(tCnt,iTrials) = NaN;
                end
            end
            SessionData.earlyDisc(SessionData.earlyDisc(:,iTrials) >= 1, iTrials) = NaN;
            SessionData.earlyDisc(SessionData.earlyDisc(:,iTrials) == 0, iTrials) = NaN;
            SessionData.lateDisc(SessionData.lateDisc(:,iTrials) >= 1, iTrials) = NaN;
            SessionData.lateDisc(SessionData.lateDisc(:,iTrials) == 0, iTrials) = NaN;
        end
        
        %% compute discrimination performance
        

        %% compute reaction times
        for iTrials = 1:SessionData.nTrials
            if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
                lLick = SessionData.RawEvents.Trial{iTrials}.Events.Port1In; %left licks
            else
                lLick = NaN; %no licks
            end
            if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Port3In') %check for licks
                rLick = SessionData.RawEvents.Trial{iTrials}.Events.Port3In; %right licks
            else
                rLick = NaN; %no licks
            end
            
            if (any(isnan(lLick)) && any(isnan(rLick))) || any(any(isnan(SessionData.RawEvents.Trial{iTrials}.States.WaitForResponse)))
                reactionTime(iTrials) = NaN; %no response
            else
                window = SessionData.RawEvents.Trial{iTrials}.States.WaitForResponse;
                lLick = lLick-window(1); lLick(lLick < 0 | lLick > (window(2) - window(1))) = [];lLick(isnan(lLick)) = inf; %compute left response during decision window
                rLick = rLick-window(1); rLick(rLick < 0 | rLick > (window(2) - window(1))) = [];rLick(isnan(rLick)) = inf; %compute right response during decision window
                
                if min(lLick) <  min(rLick)  % if left licks are first
                    reactionTime(iTrials) = lLick(1);
                elseif  min(lLick) >  min(rLick)   % if right licks are first
                    reactionTime(iTrials) = rLick(1);
                else
                    reactionTime(iTrials) = NaN;
                end
                clear window
            end
            clear  lLick rLick
        end
        
        ind = logical(SessionData.Assisted); %index for trials that were self-performed
        Performance.rTime(1,Cnt).All = nanmean(reactionTime); %mean reaction time for session
        Performance.rTime(1,Cnt).Correct = nanmean(reactionTime(SessionData.Rewarded(ind))); %mean reaction time in correct trials
        Performance.rTime(1,Cnt).False = nanmean(reactionTime(SessionData.Punished(ind))); %mean reaction time in correct trials
        
        %% combine into one larger array
        SessionData.SessionNr = repmat(Cnt,1,SessionData.nTrials); %tag all trials in current dataset with session nr
        SessionData.rTime = reactionTime; clear reactionTime
        bhv = appendBehavior(bhv,SessionData); %append into larger array
        
    end
end
disp(['Current subject: ' Animal '; Using ' num2str(Cnt) '/' num2str(size(Files,1)) ' files']);

%% check for last and high performance sessions
if highDetection > 0
    lowInd = Performance.Detection < highDetection; %find sessions with low detection performance
else
    lowInd = false(1,Cnt);
end
disp(['Rejected ' num2str(sum(lowInd)) '/' num2str(length(lowInd))  ' files for detection performance below ' num2str(highDetection*100) ' %.']);

if Cnt-sum(lowInd) < lSessions
    lSessions = Cnt-sum(lowInd);
end
disp(['Selected ' num2str(lSessions) '/' num2str(Cnt) ' remaining sessions.']);

sessionSelect = 1:lSessions; %only use last 'lSessions'
sessionSelect(lowInd(1:lSessions)) = []; %don't use low performance sessions of selected

%only use last 'lSessions' days
bhv = selectBehavior(bhv,sessionSelect); %only use trials from selecteded sessions
Performance = selectBehavior(Performance,sessionSelect); %only use performance from selecteded sessions

DayNr = DayNr(sessionSelect) - DayNr(end) + 1; %convert Days into relative values, starting at 1 for first training day
SessionNr = length(DayNr):-1:1;
fDay = Performance.Date{end}; %first date in dataset
lDay = Performance.Date{1}; %last date in dataset
disp(['First date: ' fDay]);
disp(['Last date: ' lDay]);

%% compute performance for discrimination cases (using DistStim)
[~,distStim] = histcounts((bhv.allDisc(~isnan(bhv.allDisc))),10);
distStim(end) = [];

aInd = [true false false];
vInd = [false true false];
avInd = [true true false];
    
for iDist = 1:length(distStim)
    
    ind1 = mean(bhv.allDisc(avInd,:),1) >= distStim(iDist) & mean(bhv.allDisc(avInd,:),1) < (distStim(iDist) + mean(diff(distStim))) & bhv.Assisted;
    ind2 = bhv.allDisc(vInd,:) >= distStim(iDist) & bhv.allDisc(vInd,:) < (distStim(iDist) + mean(diff(distStim))) & bhv.Assisted;
    ind3 = bhv.allDisc(aInd,:) >= distStim(iDist) & bhv.allDisc(aInd,:) < (distStim(iDist) + mean(diff(distStim))) & bhv.Assisted;
    allInd = ind1 | ind2 | ind3;
    
    [discPerf{1}(iDist,1),distConv{1}(:,iDist,1),tCount{1}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,allInd); %compute performance and error for all trials
    [discPerf{1}(iDist,2),distConv{1}(:,iDist,2),tCount{1}(iDist,2)] = computeBehavior(bhv,maxTrialCnt,allInd & bhv.CorrectSide == 1); %compute performance and error for left trials
    [discPerf{1}(iDist,3),distConv{1}(:,iDist,3),tCount{1}(iDist,3)] = computeBehavior(bhv,maxTrialCnt,allInd & bhv.CorrectSide == 2); %compute performance and error for right trials
    
    % isolate different modalites
    [modDiscPerf{1}(iDist,3),modDistConv{1}(:,iDist,3),modtCount{1}(iDist,3)] = computeBehavior(bhv,maxTrialCnt,ind1); %mixed trials
    [modDiscPerf{1}(iDist,1),modDistConv{1}(:,iDist,1),modtCount{1}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind2); %vision trials
    [modDiscPerf{1}(iDist,2),modDistConv{1}(:,iDist,2),modtCount{1}(iDist,2)] = computeBehavior(bhv,maxTrialCnt,ind3); %audio trials

end

%% Overview figure for sessions
figure('name',[Animal ' - Learning curves; Start date: ' fDay ' ; End date: ' lDay])
subplot(2,2,1)
Data = Performance.Detection;
plot(SessionNr,Data,'-ok','linewidth',2);hold on
plot(SessionNr(~isnan(Data)),smooth(Data(~isnan(Data))),'--k','linewidth',2)

line([0 length(SessionNr)+1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
title([Animal ' - Detection Performance']); xlabel('#sessions','fontsize',15); ylabel('performance (%)','fontsize',15);
axis square; set(gca,'FontSize',12); ylim([0.4 1.05]);
% legend({'Combined'},'Location','SouthEast')

subplot(2,2,2)
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
subplot(2,2,3); hold on
line([-0.2 1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])

plotDist = distStim(~isnan(discPerf{1}(:,1)));
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
    
    title([Animal ' - Discrimination Performance']); xlabel('distractor ratio','fontsize',15);
    ylabel('performance (%)','fontsize',15);
    axis square; ylim([0.35 1]); xlim([-0.1 1.1]); set(gca,'FontSize',12);
    
    subplot(2,2,4); hold on
    line([-0.2 1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    line([0.5 0.5],[-0.05 1.05],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    
    lrRatio = [plotDist ones(1,length(plotDist))]./[plotDist+1 fliplr(plotDist)+1]; %amount of right pulses as percentage of absolute pulsecount
    leftPerf = [plotPerf(:,2);flipud(1-plotPerf(:,3))]';
    leftDistConv = [distConv{1}(:,~isnan(discPerf{1}(:,1)),2) 1-fliplr(distConv{1}(:,~isnan(discPerf{1}(:,1)),3))];
    leftError(1,:) = leftPerf - leftDistConv(1,:);
    leftError(2,:) = leftDistConv(2,:) - leftPerf;
    
    errorbar(lrRatio,leftPerf,leftError(1,:),leftError(2,:), ...
        '-ok','Markersize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2);
    
    title([Animal ' - Discrimination Performance']); xlabel('right/left ratio','fontsize',15); ylabel('going left (%)','fontsize',15);
    axis square; ylim([-0.05 1.05]); xlim([-0.1 1.1]); set(gca,'FontSize',12);
end

%% compare performance at different modalities
if ~isempty(modDiscPerf{1}(~isnan(modDiscPerf{1}(:,1:2)))) %single modality trials in dataset
    
    figure;hold on
    plotError = modDistConv{1};
    plotDist = distStim;
    plotPerf = modDiscPerf{1};
    
%     plotError = modDistConv{1}(:,~isnan(modDiscPerf{1}));
%     plotDist = distStim(~isnan(modDiscPerf{1}));
%     plotPerf = modDiscPerf{1}(~isnan(modDiscPerf{1}));
    
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
    
    title([Animal ' - Single modalities']); xlabel('distractor ratio','fontsize',15); ylabel('performance (%)','fontsize',15);
    axis square; ylim([0.35 1.05]); xlim([-0.1 1.1]); set(gca,'FontSize',12);
    legend({'Vision' 'Audio' 'Multisensory'})
    line([-0.2 1.1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    line([0.5 0.5],[-0.05 1.1],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    
end