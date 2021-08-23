function [bhv] = SpatialDisc_DreaddTest(Animal,path,salineID,pharmaID,splitTrials,lSessions)
% Analze behavioral data from SpatialDisc paradigm to get basic readout of animal performance. 
% Reads all available files in 'path', ignoring data files that contain
% less then 100 trials. Saline ID is the trial ID that has to be set to
% identify saline controls (default is 2). pharmaID is the trial ID for CNO trials (default is 3).
% 
% Optional inputs: 
% lSessions: Only use last 'lSessions' sessions for behavioral analysis.
% splitTrials: Only use first and last 'splitTrials' trials and analyse
%              separately to check if performance changed over the course of a session.

%% check optional input
if ~exist('salineID','var') || isempty(salineID)
    salineID = 2;
end

if ~exist('pharmaID','var') || isempty(pharmaID)
    pharmaID = 3;
end

if ~exist('splitTrials','var') || isempty(splitTrials)
    pharmaID = 200;
end

if ~exist('lSessions','var') || isempty(lSessions)
    lSessions = inf;
end

%% assign some variables
paradigm = 'SpatialDisc';
Files = ls([path '\' Animal '\' paradigm '\Session Data\*' Animal '*']); %behavioral files in correct path
path = [path Animal '\' paradigm '\Session Data\']; %folder with behavioral data
maxTrialCnt = 500; %maximum trials per datapoint

cDate = datenum(Files(:,length([Animal '_' paradigm '_'])+1:length([Animal '_' paradigm '_'])+10)); %isolate dates from Filenames
[cDate,ind] = sort(cDate,'descend'); %sort counts to get the order of files to days correct. Newest file should be first in the list.
Files = Files(ind,:); %adjust order of filenames to get it to be chronological

%% load data
bhv = []; Cnt = 0; 
for iFiles = 1:size(Files,1)
    
    load([path Files(iFiles,:)]); %load current bhv file
    useData = isfield(SessionData,'nTrials') && SessionData.nTrials > 100 && any(SessionData.MarkerCodes == salineID | SessionData.MarkerCodes == pharmaID); % if file contains at least 100 trials
    
    if useData
        Cnt = Cnt+1;
        DayNr(Cnt) = cDate(iFiles);
        SessionData.MarkerCodes(1:end) = SessionData.MarkerCodes(1); %set marker code to all trials, so they can be separated later
        SessionData.TrialID = zeros(1,SessionData.nTrials);
        
        if SessionData.nTrials > splitTrials
            SessionData.TrialID(1:splitTrials) = 1; %give ID to first splitTrials trials to analyse performance separately.
            SessionData.TrialID(end-splitTrials+1:end) = SessionData.TrialID(end-splitTrials+1:end) + 2; %give ID to last splitTrials trials to analyse performance separately.
            SessionData.TrialID(SessionData.TrialID > 2) = 1; %make sure the first splitTrials trials are tagged correctly.
        end
        
        % combine into one larger array
        SessionData.SessionNr = repmat(Cnt,1,SessionData.nTrials); %tag all trials in current dataset with session nr
        bhv = appendBehavior(bhv,SessionData); %append into larger array
    end
end
disp(['Current subject: ' Animal '; Using ' num2str(Cnt) '/' num2str(size(Files,1)) ' files']);

%% check for last sessions
if Cnt < lSessions
    lSessions = Cnt;
end
sessionSelect = 1:lSessions; %only use last 'lSessions'
bhv = selectBehavior(bhv,sessionSelect); %only use trials from selecteded sessions
disp(['Selected ' num2str(lSessions) '/' num2str(Cnt) ' remaining sessions.']);

DayNr = DayNr(sessionSelect) - DayNr(end) + 1; %convert Days into relative values, starting at 1 for first training day
SessionNr = length(DayNr):-1:1;

%% compute performance for discrimination cases during saline or pharma injections
distFractions = bhv.DistStim ./ bhv.TargStim;
distStim = [unique(distFractions) fliplr(unique(distFractions))];
Cnt = 1;

for iDist = 1:length(distStim)
    
    if iDist == length(distStim)/2 +1
        Cnt = Cnt+1;
    end
    
    ind = distFractions == distStim(iDist) & bhv.CorrectSide == Cnt & bhv.MarkerCodes == salineID; %saline trials
    [salinePerf{1}(iDist),salineDist{1}(:,iDist,1),salineCnt{1}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind); %compute performance and error for all trials
    
    ind = distFractions == distStim(iDist) & bhv.CorrectSide == Cnt  & bhv.MarkerCodes == salineID & bhv.TrialID == 1; %saline trials, early trials
    [salinePerf{2}(iDist),salineDist{2}(:,iDist,1),salineCnt{2}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind); %compute performance and error for all trials
    
    ind = distFractions == distStim(iDist) & bhv.CorrectSide == Cnt  & bhv.MarkerCodes == salineID & bhv.TrialID == 2; %saline trials, early trials
    [salinePerf{3}(iDist),salineDist{3}(:,iDist,1),salineCnt{3}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind); %compute performance and error for all trials
    
    ind = distFractions == distStim(iDist) & bhv.CorrectSide == Cnt  & bhv.MarkerCodes == pharmaID; %pharma trials
    [pharmaPerf{1}(iDist),pharmaDist{1}(:,iDist,1),pharmaCnt{1}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind); %compute performance and error for all trials
    
    ind = distFractions == distStim(iDist) & bhv.CorrectSide == Cnt  & bhv.MarkerCodes == pharmaID & bhv.TrialID == 1; %pharma trials, early trials
    [pharmaPerf{2}(iDist),pharmaDist{2}(:,iDist,1),pharmaCnt{2}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind); %compute performance and error for all trials
    
    ind = distFractions == distStim(iDist) & bhv.CorrectSide == Cnt  & bhv.MarkerCodes == pharmaID & bhv.TrialID == 2; %pharma trials, early trials
    [pharmaPerf{3}(iDist),pharmaDist{3}(:,iDist,1),pharmaCnt{3}(iDist,1)] = computeBehavior(bhv,maxTrialCnt,ind); %compute performance and error for all trials
    
end

for x = 1:3
    salinePerf{x}(1:length(distStim)/2) = 1 - salinePerf{x}(1:length(distStim)/2);
    salineDist{x}(:,1:length(distStim)/2) = 1 - salineDist{x}(:,1:length(distStim)/2);
    
    pharmaPerf{x}(1:length(distStim)/2) = 1 - pharmaPerf{x}(1:length(distStim)/2);
    pharmaDist{x}(:,1:length(distStim)/2) = 1 - pharmaDist{x}(:,1:length(distStim)/2);
end

%% overview for distractor performance
plotDist = distStim(1:length(distStim)/2);
plotDist = [plotDist ones(1,length(plotDist))]./[plotDist+1 fliplr(plotDist)+1]; %amount of right pulses as percentage of absolute pulsecount
titelText = {'Total' 'Early' 'Late'};

figure;
for x = 1:3
    subplot(1,3,x); hold on
    line([-0.1 1.1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    line([0.5 0.5],[-0.05 1.05],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    
    Est = [pharmaPerf{x}(end) 0.5  0.1];
    f = @(p,x) p(1) ./ (1 + exp(((1-x)-p(2))/p(3)));
    wFit1 = fitnlm(plotDist,(pharmaPerf{x}),f,Est,'Weight',(pharmaCnt{1}.*~(pharmaCnt{1}<0.45))'+1);
    wFit2 = fitnlm(plotDist,(salinePerf{x}),f,Est,'Weight',(pharmaCnt{1}.*~(pharmaCnt{1}<0.45))'+1);
  
    
    errorbar(plotDist,pharmaPerf{x},pharmaPerf{x} - pharmaDist{x}(1,:),pharmaDist{x}(2,:) - pharmaPerf{x}, ...
        'ok','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2)
    line(xx,predict(wFit1,xx),'Linewidth',2,'color','k'); % plot fit
    
    errorbar(plotDist,salinePerf{x},salinePerf{x} - salineDist{x}(1,:),salineDist{x}(2,:) - salinePerf{x}, ...
        'og','Markersize',5,'MarkerEdgeColor','g','MarkerFaceColor','w','linewidth',2)
    line(xx,predict(wFit2,xx),'Linewidth',2,'color','g'); % plot fit

    title([Animal ' - ' titelText{x} ' Discrimination']); xlabel('right/left ratio','fontsize',15); ylabel('going right (%)','fontsize',15);
    axis square; ylim([-0.05 1.05]); xlim([-0.1 1.1]); set(gca,'FontSize',12);
    
end
end