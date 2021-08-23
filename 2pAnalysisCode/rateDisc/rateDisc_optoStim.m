function [Performance,bhv] = rateDisc_optoStim(Animal,cPath,lSessions,highDetection,discOnly)
% Analyze behavioral data from rate discrimination task to test for the
% impact of optogenetic manipulation.
%
% Optional inputs:
% lSessions: Only use sessions after 'lSessions' which should be datestring that can be converted into a datenumber.
% highDetection: Only use sessions at which detection was at 90% or higher.
% discOnly: Only use sessions at which distractor stimuli were presented.

%% check optional input
if ~exist('lSessions','var') || isempty(lSessions)
    lSessions = inf;
end

if ~exist('highDetection','var') || isempty(highDetection)
    highDetection = 0;
end

if ~exist('discOnly','var') || isempty(discOnly)
    discOnly = false;
end

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end
cropLength = 50; %nr of trials in the beginning/end of session that get cropped if performance is below 60%.

%% get files and date for each recording
paradigm = 'SpatialDisc';
for iChecks = 1:10 %check for files repeatedly. Sometimes the server needs a moment to be indexed correctly
    Files = ls([cPath '\' Animal '\' paradigm '\Session Data\*' Animal '*']); %behavioral files in correct cPath
    if ~isempty(Files)
        break;
    end
    pause(0.1);
end

cPath = [cPath Animal '\' paradigm '\Session Data\']; %folder with behavioral data

cDate = datenum(Files(:,length([Animal '_' paradigm '_'])+1:length([Animal '_' paradigm '_'])+10)); %isolate dates from Filenames
cDate = cDate + (str2num(Files(:,length([Animal '_' paradigm '_'])+19:end-4))*0.01); %add session nr to timestamp

[cDate,ind] = sort(cDate,'descend'); %sort counts to get the order of files to days correct. Newest file should be first in the list.
cDate = floor(cDate);
Files = Files(ind,:); %adjust order of filenames to get it to be chronological

if any(~isinf(lSessions))
    cIdx = cDate >= datenum(lSessions); %reject session before requested date
    cDate = cDate(cIdx);
    Files = Files(cIdx,:);
end

%% load data
bhv = []; Performance = []; Cnt = 0;

for iFiles = 1:size(Files,1)
    try
        clear SessionData
        load([cPath Files(iFiles,:)], 'SessionData'); %load current bhv file
        if isfield(SessionData,'Rewarded')
            SessionData.Rewarded = logical(SessionData.Rewarded);
        end
        
        useData = false;
        if isfield(SessionData, 'optoDur')
            if sum(SessionData.optoDur > 0) > 0
                useData = length(SessionData.Rewarded) > 100; % if file contains at least 100 trials;
                
                if discOnly && useData %only discrimination sessions
                    useData = sum(SessionData.DistStim > 0) > 100; % if file contains at least 100 trials and different distractor rates;
                end
            end
            
            if isempty(SessionData.Notes{1}) && useData
                disp(['No fiber location in recording: ' Files(iFiles,:)]);
                useData = false;
            end
        end
    catch
        useData = false;
    end
    
    if useData
        if isfield(SessionData.TrialSettings(end), 'optoPower')
            SessionData.optoPower = [SessionData.TrialSettings.optoPower];
        else
            SessionData.optoPower = NaN(1,length(SessionData.Rewarded));
        end
    end

    if useData
        Cnt = Cnt+1;
        DayNr(Cnt) = cDate(iFiles);
        Performance.Animal{1,Cnt} = Animal;
        Performance.Date{1,Cnt} = datestr(cDate(iFiles));
        if ischar(SessionData.Notes{1})
            Performance.Notes{1,Cnt} = strtrim(SessionData.Notes{1}(1,:));
        elseif ~isempty(SessionData.Notes{1})
            Performance.Notes{1,Cnt} = strtrim(SessionData.Notes{1}{1});
        end
        if strcmpi(Performance.Notes{1,Cnt}, 'frontal')
            SessionData.stimLocation = ones(1, length(SessionData.Rewarded)) * 1;
        elseif strcmpi(Performance.Notes{1,Cnt}, 'parietal')
            SessionData.stimLocation = ones(1, length(SessionData.Rewarded)) * 2;
        elseif strcmpi(Performance.Notes{1,Cnt}, 'v1')
            SessionData.stimLocation = ones(1, length(SessionData.Rewarded)) * 3;
        else
            SessionData.stimLocation = zeros(1, length(SessionData.Rewarded));
        end
        
        SessionData.date = repmat(cDate(iFiles), 1, length(SessionData.Rewarded));
        
        % cut beginning/end of session is performance is too low
        if mean(SessionData.Rewarded(1:cropLength)) < 0.6 %don't use intiial trials if performance is too low
            SessionData = selectBehaviorTrials(SessionData, cropLength + 1 : size(SessionData.Rewarded,2));
        end
        if mean(SessionData.Rewarded(end-cropLength:end)) < 0.6 %don't use last trials if performance is too low
            SessionData = selectBehaviorTrials(SessionData, 1 : size(SessionData.Rewarded,2) - cropLength);
        end
        
        ind = ~SessionData.DidNotChoose & ~SessionData.DidNotLever & logical(SessionData.Assisted) & SessionData.DistStim == 0; %only use active detection trials
        dInd = ind & SessionData.optoDur == 0;
        Performance.Detection(Cnt) = sum(SessionData.Rewarded(dInd))/sum(SessionData.Rewarded(dInd)+SessionData.Punished(dInd)); %peformance for self-performed trials
        
        optCnt = 0;
        for optoSide = unique(SessionData.optoSide(~isnan(SessionData.optoSide)))
            for optoTime = unique(SessionData.optoType(~isnan(SessionData.optoType)))
                optCnt = optCnt + 1;
                cInd = ind & SessionData.optoDur > 0 & SessionData.optoSide == optoSide & SessionData.optoType == optoTime;
                Performance.optoDetection{Cnt}(optCnt) = sum(SessionData.Rewarded(cInd))/sum(SessionData.Rewarded(cInd)+SessionData.Punished(cInd)); %peformance for self-performed trials
                Performance.optoSide{Cnt}(optCnt) = optoSide; %current side
                Performance.optoTime{Cnt}(optCnt) = optoTime; %current time
            end
        end
        
        %% combine into one larger array
        SessionData.SessionNr = repmat(Cnt,1,SessionData.nTrials); %tag all trials in current dataset with session nr
        bhv = appendBehavior(bhv,SessionData); %append into larger array
        
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
    
    fDay = Performance.Date{end}; %first date in dataset
    lDay = Performance.Date{1}; %last date in dataset
    disp(['First date: ' fDay]);
    disp(['Last date: ' lDay]);
end


% old code
%
%     %% compute performance for all data combined
%     for iOpto = 1:2
%         for iMod = 1 : 4
%             if iMod == 4
%                 ind = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted); %only use active trials, all modalities
%             else
%                 ind = ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.StimType == modId(iMod); %only use active trials
%             end
%             
%             if iOpto == 2
%                 ind = ind & bhv.optoDur > 0 & bhv.optoSide == optoSide & bhv.optoType == optoType;
%             end
%             
%             allPerf.SelfPerformed(iOpto, iMod) = sum(bhv.Rewarded(ind))/sum(bhv.Rewarded(ind)+bhv.Punished(ind)); %peformance for self-performed trials
%             dInd = logical(bhv.Assisted & bhv.DistStim == 0); %index for trials that were detection only
%             allPerf.Detection(iOpto, iMod) = sum(bhv.Rewarded(ind & dInd))/sum(bhv.Rewarded(ind & dInd)+bhv.Punished(ind & dInd)); %peformance for detection trials
%             allPerf.Discrimination(iOpto, iMod) = sum(bhv.Rewarded(ind & ~dInd))/sum(bhv.Rewarded(ind & ~dInd)+bhv.Punished(ind & ~dInd)); %peformance for discrimination trials
%             
%             lInd = bhv.CorrectSide == 1; %index for left-choice trials
%             allPerf.LeftPerformed(iOpto, iMod) = sum(bhv.Rewarded(ind & lInd))/sum(bhv.Rewarded(ind & lInd)+bhv.Punished(ind & lInd));
%             allPerf.RightPerformed(iOpto, iMod) = sum(bhv.Rewarded(ind & ~lInd))/sum(bhv.Rewarded(ind & ~lInd)+bhv.Punished(ind & ~lInd));
%             
%             %% compute discrimination performance and stats
%             if ~isnan(allPerf.Discrimination(iOpto, iMod))
%                 rInd = (bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded); %right-choice trials
%                 rInd = rInd(ind);
%                 
%                 tCnt = 0;
%                 eventCnt = zeros(2,sum(ind));
%                 for iTrials = find(ind)
%                     tCnt = tCnt +1;
%                     [left, right] = Behavior_getStimEvent(bhv.StimType(iTrials), bhv.stimEvents{iTrials});
%                     eventCnt(1,tCnt) = length(right);
%                     eventCnt(2,tCnt) = length(left) + length(right);
%                 end
%                 
%                 [nTrials, distRatio, trialInd]  = histcounts(eventCnt(1,:) ./ eventCnt(2,:), binSize); %get ratio between rightward and absolute number of events
%                 distRatio = distRatio + diff(distRatio(1:2))/2; distRatio(end) = [];
%                 for iBins = 1:length(nTrials)
%                     rightChoice(iBins) = sum(rInd(trialInd == iBins)); %get number of rightward trials for each difficulty
%                 end
%                 
%                 [params, h1, ~, cFit] = Behavior_fitPalamedes(distRatio, rightChoice, nTrials, showPlot, true); %complete discrimination parameters
%                 
%                 if ~isempty(h1)
%                     if iOpto == 1
%                         title(h1.data.Parent,[Animal ' - ' modLabels{iMod} ' - All non-stimulated trials'])
%                     elseif iOpto == 2
%                         title(h1.data.Parent,[Animal ' - ' modLabels{iMod} ' - All optogenetic trials'])
%                     end
%                     ylim(h1.data.Parent, [0 1]);
%                 end
%                 
%                 allPerf.cFit{iOpto, iMod} = cFit;
%                 
%                 fNames = fieldnames(params);
%                 for iFields = 1 : length(fNames)
%                     allPerf.(fNames{iFields})(iMod) = params.(fNames{iFields});
%                 end
%                 
%                 % discrimination only
%                 [params, ~, ~, cFit] = Behavior_fitPalamedes(distRatio(2:end-1), rightChoice(2:end-1), nTrials(2:end-1)); %complete discrimination parameters
%                 allPerf.disc_cFit{iOpto, iMod} = cFit;
%                 if showPlot
%                     hold(h1.fit.Parent, 'on');
%                     h2 = plot(h1.fit.Parent,cFit(1,:),cFit(2,:),'color','r', 'linewidth', 4);
%                     uistack(h2,'bottom');
%                     hold(h1.fit.Parent, 'off');
%                 end
%                 
%                 fNames = fieldnames(params);
%                 for iFields = 1 : length(fNames)
%                     allPerf.(['disc_' fNames{iFields}])(iMod) = params.(fNames{iFields});
%                 end
%                 
%                 allDisc{iOpto}(iMod,:) = rightChoice ./ nTrials;
%                 [upper, lower] = Behavior_wilsonError(rightChoice, nTrials);
%                 allDiscUpper{iOpto}(iMod,:) = upper;
%                 allDiscLower{iOpto}(iMod,:) = lower;
%             end
%         end
%     end
%     
%     %% make figure
%     if showPlot
%         figure;
%         for iMod = 1 : 4
%             subplot(2,2,iMod);
%             plot(allPerf.cFit{1, iMod}(1,:), allPerf.cFit{1, iMod}(2,:),'color','k', 'linewidth', 4); hold on;
%             plot(allPerf.cFit{2, iMod}(1,:), allPerf.cFit{2, iMod}(2,:),'color','b', 'linewidth', 4);
%             
%             errorbar(distRatio, allDisc{1}(iMod,:), allDisc{1}(iMod,:) - allDiscLower{1}(iMod,:), allDiscUpper{1}(iMod,:) - allDisc{1}(iMod,:),'ok', 'MarkerFaceColor','w','linewidth',2);
%             errorbar(distRatio, allDisc{2}(iMod,:), allDisc{2}(iMod,:) - allDiscLower{2}(iMod,:), allDiscUpper{2}(iMod,:) - allDisc{2}(iMod,:),'ob', 'MarkerFaceColor','w','linewidth',2);
%             
%             ylim([0, 1]); xlim([0, 1]); axis square
%             xlabel('pulseRatio'); ylabel('ChoseRight');
%             title(modLabels{iMod})
%         end
%     end