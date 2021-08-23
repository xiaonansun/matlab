function rateDisc_runDecoderSegs(animal)
%% some basic variables
trialsPerBin = 200;
decTypes = {'all','correct','error','choice','stim'};
cMods = [0 2 4 6];

%% path and data
if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for global dimensions

load([bPath 'trialInfo.mat'],'bhvTrials','recs','trialCnt');
load([bPath animal '_blockBhv.mat'],'bhv');

allTrials = [0 cumsum(trialCnt)]; %trials per sessions

% for each window, find the required amount of cases when combining choice, stimulus and correct/error
leftIdx = (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left (choice)
corrIdx = bhv.Rewarded; %rewarded trials

%% loop through conditions and send jobs to hpc system
for iMod = 1 : length(cMods)
    if cMods(iMod) == 0
        modIdx = true(1,length(bhv.StimType)); %use all trials
    else
        modIdx =  bhv.StimType == cMods(iMod); %select modality
    end
    
    for iDec = 1 : length(decTypes)
        % get combined cases to determine window size
        binFac = 1.5; %by default, use three times more trials as required to get trials in the delay
        clear tCombs
        if strcmpi(decTypes{iDec}, 'all')
            tCombs{1} = find(leftIdx & modIdx);
            tCombs{2} = find(~leftIdx & modIdx);
        elseif strcmpi(decTypes{iDec}, 'correct')
            tCombs{1} = find(leftIdx & modIdx & corrIdx);
            tCombs{2} = find(~leftIdx & modIdx & corrIdx);
        elseif strcmpi(decTypes{iDec}, 'error')
            tCombs{1} = find(leftIdx & modIdx & ~corrIdx);
            tCombs{2} = find(~leftIdx & modIdx & ~corrIdx);
        elseif strcmpi(decTypes{iDec}, 'choice') || strcmpi(decTypes{iDec}, 'stim')
            tCombs{1} = find(corrIdx & leftIdx & modIdx);
            tCombs{2} = find(corrIdx & ~leftIdx & modIdx);
            tCombs{3} = find(~corrIdx & ~leftIdx & modIdx);
            tCombs{4} = find(~corrIdx & leftIdx & modIdx);
            binFac = binFac/2; %this needs to be adjusted to account for using 4 cases instead of 2
        end
        
        % determine total number of steps
        lastTrial = 0; breaker = false;
        for iSteps = 1 : 10000
            cIdx = [lastTrial lastTrial];
            for x = 1 : length(tCombs)
                temp = tCombs{x}(tCombs{x} > lastTrial); %find the case that needs the most trials to reach the required case count
                if length(temp) > ceil(trialsPerBin*binFac)
                    y = temp(ceil(trialsPerBin*binFac));
                elseif ~isempty(temp)
                    y = temp(end);
                    breaker = true;
                end
                if y > cIdx(2)
                    cIdx(2) = y; %use selected case to find the end of the trialindex
                end
            end
            lastTrial = lastTrial + round(diff(cIdx) / 5); %move with 80% overlap
            if breaker; break; end
        end
        stepCnt = iSteps; %total number of steps
        
        %% run through steps and send to system
        for iSteps = 1 : stepCnt
            disp([animal ' - ' decTypes{iDec} ' - iMod: ' num2str(cMods(iMod)) ' - Step ' num2str(iSteps) '/' num2str(stepCnt)]);
            cLine = ['qsub -l m_mem_free=3G -pe threads 8 -binding linear:8 choiceDecode.sh ' ...
                animal ' ' num2str(iSteps) ' ' decTypes{iDec} ' ' num2str(cMods(iMod)) ' ' num2str(trialsPerBin)];
            system(cLine);
        end
    end
end
