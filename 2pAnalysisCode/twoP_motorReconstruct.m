function [recV, semV, allBeta, recLabels, dataPath, allModIdx, allSideIdx, alignIdx, allRecIdx, baseLength, frames, stimTimes, cellCnt] = twoP_motorReconstruct(cMod)
% code to create reconstructed 2pimaging data based on motor regressors. cMod
% informs over which animals to use ('Visual', 'Audio' or 'All'). Last
% dimension of recV is different modalities, 1 is corr. vision, 2 is corr. audio, 
% 3 is all corr. trials. 4:6 is the same thing but for all trials.

%% select data sets
dataOverview = twoP_delayDecRecordings; %get overview for all recordings
[~, motorLabels] = delayDecRecordings; %get index for motor regressors
opMotorLabels = {'lLick' 'rLick' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'}; %operant motor regressors

animals = dataOverview(:,1);
Cnt = 0;
baseLength = inf;
postLength = inf;

for iAnimals = 1:length(animals)
    
    fPath = dataOverview{iAnimals,5};
    if strcmpi(dataOverview{iAnimals,2},cMod) || strcmpi(cMod,'all') && exist([fPath 'dimBeta.mat'],'file')    
        %% load data
        Cnt = Cnt +1;
        fPath = dataOverview{iAnimals,5};
        dataPath{Cnt} = fPath; %store current data path
        load([fPath 'dimBeta.mat'],'dimBeta');
        load([fPath 'regData.mat'],'fullR','trialIdx','recIdx','idx','recLabels'); %load model data
        cellCnt(Cnt+1) = size(dimBeta,2);
        allBeta{Cnt} = dimBeta;
        allRecIdx{Cnt} = recIdx(~idx);
        
        %% get frametimes to determine #frames / trial and construct modality indices
        load([fPath 'data.mat'],'data');
        load([fPath 'interpVc.mat'],'frames');
        trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
        
        % load behavior and get modality indices
        bhvFile = strsplit(fPath,filesep);
        bhvFile = dir([fPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
        load([fPath bhvFile.name]);
        
        bTrials = data.trialNumbers;
        bTrials(~ismember(data.trialNumbers,data.bhvTrials)) = []; %don't use trials that have problems with trial onset times
        bTrials(SessionData.DidNotChoose(bTrials) | SessionData.DidNotLever(bTrials) | ~SessionData.Assisted(bTrials)) = []; %don't use unperformed/assisted trials
        
        for iTrials = 1 : length(SessionData.Rewarded)
            try
                leverTimes = [reshape(SessionData.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
                    reshape(SessionData.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
                    reshape(SessionData.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
                stimGrab = leverTimes(find(leverTimes == SessionData.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
                stimTimes{Cnt}(iTrials) = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab; %time of stimulus onset - measured from soundcard
            catch
                stimTimes{Cnt}(iTrials) = NaN;
            end
        end
        clear stimGrab leverTimes
        
        %% realign data so baseline is aligned to handle grab and poststim to stimulus
        stimOn = sum(fullR(:,ismember(recIdx(~idx),find(ismember(recLabels,{'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'})))),2); %index for stimulus onset in all trials
        stimOn = find([0;diff(stimOn)] > 0.5) - 1;

        % index for baseline (time before first possible stimulus onset)
        baseLength = min([min(unique(rem(stimOn,frames)))-1 baseLength]);
        baseIdx = repmat((0:frames:size(fullR,1)-1)',1,baseLength);
        baseIdx = bsxfun(@plus,baseIdx, 1 : baseLength);
        baseIdx = baseIdx(:);

        % index for post stimulus time
        postLength = min([frames - max(unique(rem(stimOn,frames))) postLength]); %shortest possible poststim duration
        stimIdx = repmat(stimOn, 1, postLength);
        stimIdx = bsxfun(@plus,stimIdx, 0 : postLength-1);
        stimIdx = stimIdx(:);
        
        alignIdx{Cnt} = sort([baseIdx;stimIdx]);
        fullR = bsxfun(@minus, fullR, mean(fullR, 1)); %make sure design matrix is zero-mean
        fullR = fullR(alignIdx{Cnt},:); %reduce fullR to only include aligned baseline and poststim data
        frames = postLength + baseLength; %new single trial duration in frames
        
        %% build lots of indices
        sucInd = SessionData.Rewarded(bTrials) & SessionData.Assisted(bTrials); %find succesful unisensory trials
        modIdx(1,:) = reshape(repmat(SessionData.StimType(bTrials) == 1 & sucInd,frames,1),[],1); % correct visual trials
        modIdx(2,:) = reshape(repmat(SessionData.StimType(bTrials) == 2 & sucInd,frames,1),[],1); % correct audio trials
        modIdx(3,:) = sum(modIdx(1:2,:)); % all correct trials
        modIdx(4,:) = reshape(repmat(SessionData.StimType(bTrials) == 1, frames,1),[],1); % all visual trials
        modIdx(5,:) = reshape(repmat(SessionData.StimType(bTrials) == 2, frames,1),[],1); % all audio trials        
        modIdx(6,:) = sum(modIdx(4:5,:)); % all trials
        modIdx = reshape(modIdx,size(modIdx,1),frames,[]);
        modIdx = modIdx(:,:,trialIdx); %remove non-used trials
        modIdx = reshape(modIdx,size(modIdx,1),[]);

        sideIdx(1,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 1 & sucInd,frames,1),[],1); % correct left trials
        sideIdx(2,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 2 & sucInd,frames,1),[],1); % correct right trials
        sideIdx(3,:) = sum(sideIdx(1:2,:)); % all correct trials
        sideIdx(4,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 1, frames,1),[],1); % left trials
        sideIdx(5,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 2, frames,1),[],1); % right trials
        sideIdx(6,:) = sum(sideIdx(4:5,:)); % all trials    
        sideIdx = reshape(sideIdx,size(sideIdx,1),frames,[]);
        sideIdx = sideIdx(:,:,trialIdx); %remove non-used trials
        sideIdx = reshape(sideIdx,size(sideIdx,1),[]);
        
        %% cycle through regressors and reconstruct each one        
        for iRegs = 1 : length(recLabels) + 4
            if iRegs <= length(recLabels)
                cInd = ismember(recIdx(~idx), find(ismember(recLabels,recLabels{iRegs}))); %find current regressors
            elseif iRegs == length(recLabels) + 1
                cInd = ismember(recIdx(~idx), find(~ismember(recLabels,motorLabels))); %find task regressors
            elseif iRegs == length(recLabels) + 2
                cInd = ismember(recIdx(~idx), find(ismember(recLabels,opMotorLabels))); %find operant motor regressors
            elseif iRegs == length(recLabels) + 3
                cInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels) & ~ismember(recLabels,opMotorLabels))); %find spont. motor regressors
            elseif iRegs == length(recLabels) + 4
                cInd = ismember(recIdx(~idx), find(ismember(recLabels,recLabels))); %find all regressors
            end
            
            data = fullR(:, cInd); %current regressors
            cBeta = dimBeta(cInd,:); %current weights
            Vm = (data * cBeta)'; %model Vc data
            
            if iAnimals == 1
                recV{iRegs} = [];
                semV{iRegs} = [];
            end
            
            for iMod = 1:8 %different modalities, 1 is corr. vision, 2 is corr. audio, 3 is all corr. trials. 4:6 is the same thing but for all trials. 7 is all leftward trials, 8 is all rightward trials
                
                if iMod == 7
                    temp = reshape(Vm(:,sideIdx(1,:)),size(Vm,1),frames,[]); %all correct leftward trials
                    cError = std(temp,[],3);
                elseif iMod == 8
                    temp = reshape(Vm(:,sideIdx(2,:)),size(Vm,1),frames,[]); %all correct rigthward trials
                    cError = std(temp,[],3);
                else
                    temp = reshape(Vm(:,modIdx(iMod,:)),size(Vm,1),frames,[]);
                    cError = sem(temp,3);
                end
                temp = mean(temp,3);
                if size(recV{iRegs},2) > size(temp,2)
                    recV{iRegs}(:,size(temp,2)+1:end,:) = [];
                    semV{iRegs}(:,size(cError,2)+1:end,:) = [];
                end
                if iMod == 6 && any(isnan(temp(:)))
                    error('error');
                end
                recV{iRegs}(sum(cellCnt(1:Cnt))+1 : sum(cellCnt(1:Cnt+1)),:,iMod) = temp; %reconstructed V for current motor regressor and modality
                semV{iRegs}(sum(cellCnt(1:Cnt))+1 : sum(cellCnt(1:Cnt+1)),:,iMod) = cError; %reconstructed V error for current motor regressor and modality
            end
        end
        
        for iMod = 1:6
            allModIdx{Cnt,iMod} = modIdx(iMod,:);
            allSideIdx{Cnt,iMod} = sideIdx(iMod,:);
        end
        clear sideIdx modIdx Vm
    end
end
cellCnt(1) = []; %first entry is not a recording, remove.
end
 
        