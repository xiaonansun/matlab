% function rateDisc_checkMoveModel
%code to asses how much neural variance is explained by movements variables in the linear encoding model across ssessions

% Convention in recIdx: 
% 1 = early audio learning
% 2 = 5 to 50 percentile audio detection
% 3 = 50 to 95 percentile audio detection
% 4 = remaining audio detection
% 5 = Audio discrimination
% 6 = Novice tactile discrimination
% 7 = 5 to 50 percentile tactile detection
% 8 = 50 to 95 percentile tactile detection
% 9 = remaining tactile detection
% 10 = tactile discrimination
% 11 = everything else (usually mixed discriminnation)

%% some variables
if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
    gPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
    bhvPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
    bhvPath = '/sonas-hs/churchland/nlsas/data/data/Behavior_Simon/';
end
load allenDorsalMapSM
allenMask = dorsalMaps.allenMask;

% this is to isolate movement events
opMoveLabels = {'lLick' 'rLick' 'lGrab' 'rGrab'};
spontMoveLabels = {'piezo' 'whisk' 'nose' 'body'};
moveShift = 16; %timepoint 0 in movement design matrx

% this is to isolate beta weights
[~, motorLabels] = delayDecRecordings;
allOp = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'}; %all operant motor regressors
allSpont = motorLabels(~ismember(motorLabels,allOp)); %all spontaneous motor regressors

%% get list of recordings and reference image and load data from all sessions
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
% animals = {'mSM63' 'mSM64' 'mSM65' 'mSM66'};
% animals = {'Fez7' 'Fez10'};
animals = {'Plex01' 'Plex02'};

fullMaps = cell(1, length(animals));
taskMaps = cell(1, length(animals));
motorMaps = cell(1, length(animals));

opMoves = cell(1, length(animals));
opMoveTrialVar = cell(1, length(animals));
spontMoves = cell(1, length(animals));
spontMoveTrialVar = cell(1, length(animals));

taskPower = cell(1, length(animals));
opMotorPower = cell(1, length(animals));
spontMotorPower = cell(1, length(animals));

allData = cell(1, length(animals));
allTrialVar = cell(1, length(animals));

expTrialVar = cell(1, length(animals));
expAvgVar = cell(1, length(animals));

for iAnimals = 1 : length(animals)
    
    animal = animals{iAnimals};
    aPath = [cPath animal filesep];
    
    [dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
    recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
    recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs, 'rateDisc'); %this sorts recordings by date
    fprintf('Basepath: %s; Found %d recordings\n', aPath, length(recs));
    
    recIdx{iAnimals} = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs); %for each recording, determine to which part of the training it belongs
    firstDate(iAnimals) = rateDisc_checkFirstRec(bhvPath, animal); %get time of first behavioral session
    
    % apply U to raw data compute temporal component
    Cnt = 0; useIdx = false(1,length(recs));
    for iRecs = 1 : find(recIdx{iAnimals} == 5,1,'last')
%     for iRecs = 1 : length(recs)
        tic
        % get model data
        fPath = [aPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        gfPath = [gPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        checker = false;
        
        try
            load([fPath 'fullCorr.mat'],'fullMap');
            load([fPath 'taskregData.mat'],'taskMap');
            load([fPath 'motorregData.mat'],'motorMap');
            load([fPath 'taskSpMotorregData.mat'],'taskSpMotorMap');
            load([fPath 'taskOpMotorregData.mat'],'taskOpMotorMap');
            load([fPath 'Vc.mat'],'U','bTrials','Sv','totalVar');
            load([fPath 'interpVc.mat'],'Vc','frames');
            load([fPath 'mask.mat'],'mask');
            load([fPath 'opts2.mat'],'opts');
            load([fPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'stdV');
            load([fPath 'BehaviorVideo' filesep 'FilteredPupil.mat'],'faceM');
            load([fPath 'regData.mat'],'rejIdx','regIdx','regLabels','trialIdx');
            load([fPath 'spontMotorregData.mat'],'spontMotorMap','spontMotorR','spontMotorIdx','spontMotorLabels');
            load([fPath 'opMotorregData.mat'],'opMotorMap','opMotorR','opMotorIdx','opMotorLabels');
            load([fPath 'taskBeta.mat'],'taskBeta');
            load([fPath 'spontMotorBeta.mat'],'spontMotorBeta');
            load([fPath 'opMotorBeta.mat'],'opMotorBeta');
            
            load([fPath 'interpVtask.mat'], 'Vtask'); %load model reconstruction data
            load([fPath 'interpVspontMotor.mat'], 'VspontMotor'); %load model reconstruction data
            load([fPath 'interpVopMotor.mat'], 'VopMotor'); %load model reconstruction data
            
            U = arrayShrink(U,mask,'merge');

            bhvFile = dir([fPath filesep animal '_SpatialDisc*.mat']);
            bhvFile = strtrim(bhvFile.name);
            load([fPath bhvFile]); %load behavior data
            trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
            bhv = selectBehaviorTrials(SessionData,bTrials(trialIdx)); %only use completed trials that are in the Vc dataset
            leftIdx = bhv.Rewarded & (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left and got rewarded

            if nanmean(fullMap(:)) > 0.8 || length(bTrials) < 50 || frames ~= 150; error('not good'); end
            checker = true;
        end
        if checker
            Cnt = Cnt +1; useIdx(iRecs) = true;
            fullMap = fullMap.^2; taskMap = taskMap.^2; motorMap = motorMap.^2;
            spontMotorMap = spontMotorMap.^2; opMotorMap = opMotorMap.^2;
            taskSpMotorMap = taskSpMotorMap.^2; taskOpMotorMap = taskOpMotorMap.^2;
            
            % compute trial-by-trial correlations
            Vc = reshape(Vc,size(Vc,1),frames,[]);
            Vtask = reshape(Vtask,size(Vtask,1),frames,[]);
            VopMotor = reshape(VopMotor,size(VopMotor,1),frames,[]);
            VspontMotor = reshape(VspontMotor,size(VspontMotor,1),frames,[]);
            
            % compute trial-by-trial and trial-averaged correlations
            [expTrialVar{iAnimals}(:,Cnt,1) , expAvgVar{iAnimals}(:,Cnt,1)] = rateDisc_computeTrialCorr(U, Vc, Vtask, mask, allenMask, opts);
            [expTrialVar{iAnimals}(:,Cnt,2) , expAvgVar{iAnimals}(:,Cnt,2)] = rateDisc_computeTrialCorr(U, Vc, VopMotor, mask, allenMask, opts);
            [expTrialVar{iAnimals}(:,Cnt,3) , expAvgVar{iAnimals}(:,Cnt,3)] = rateDisc_computeTrialCorr(U, Vc, VspontMotor, mask, allenMask, opts);
           
            %reconstruct time after averaging over all pixels and get average and variance
            Vc = reshape(Vc,size(Vc,1),[]);
            Vc = mean(U) * Vc;
            Vc = reshape(Vc, frames, []);
            Vc = Vc(:,leftIdx); %only use correct, left-choice trials
            allData{iAnimals}(Cnt,:,1) = nanmean(Vc,2); %mean over trials
            allData{iAnimals}(Cnt,:,2) = nanvar(Vc,[],2); %variance over trials in each time point
            allTrialVar{iAnimals}(Cnt) = var(nanmean(abs(Vc),1)); %trial by trial variability
            
            allVar{iAnimals}(Cnt) = totalVar; %total variance in the data
            allDim{iAnimals}(Cnt) = find(cumsum(Sv./sum(Sv)) > 0.99,1); %how many dimensions needed
            allMaps{iAnimals}(Cnt,:) = [mean(fullMap) mean(motorMap) mean(taskMap) mean(spontMotorMap) mean(opMotorMap) mean(taskSpMotorMap) mean(taskOpMotorMap)];
            
            allDates{iAnimals}(Cnt) = datenum(recs(iRecs).name);
            trialCnt{iAnimals}(Cnt) = length(bTrials);
            allMoveVar{iAnimals}(Cnt) = sum(stdV(1:200).^2);
            
            fullMap = alignAllenTransIm(arrayShrink(fullMap,mask,'split'), opts.transParams);
            fullMaps{iAnimals}(:, Cnt) = arrayShrink(fullMap(:, 1:size(allenMask,2)),allenMask,'merge');
            taskMap = alignAllenTransIm(arrayShrink(taskMap,mask,'split'), opts.transParams);
            taskMaps{iAnimals}(:, Cnt) = arrayShrink(taskMap(:, 1:size(allenMask,2)),allenMask,'merge');
            motorMap = alignAllenTransIm(arrayShrink(motorMap,mask,'split'), opts.transParams);
            motorMaps{iAnimals}(:, 1, Cnt) = arrayShrink(motorMap(:, 1:size(allenMask,2)),allenMask,'merge');
            spontMotorMap = alignAllenTransIm(arrayShrink(spontMotorMap,mask,'split'), opts.transParams);
            motorMaps{iAnimals}(:, 2, Cnt) = arrayShrink(spontMotorMap(:, 1:size(allenMask,2)),allenMask,'merge');
            opMotorMap = alignAllenTransIm(arrayShrink(opMotorMap,mask,'split'), opts.transParams);
            motorMaps{iAnimals}(:, 3, Cnt) = arrayShrink(opMotorMap(:, 1:size(allenMask,2)),allenMask,'merge');
            taskSpMotorMap = alignAllenTransIm(arrayShrink(taskSpMotorMap,mask,'split'), opts.transParams);
            motorMaps{iAnimals}(:, 4, Cnt) = arrayShrink(taskSpMotorMap(:, 1:size(allenMask,2)),allenMask,'merge');
            taskOpMotorMap = alignAllenTransIm(arrayShrink(taskOpMotorMap,mask,'split'), opts.transParams);
            motorMaps{iAnimals}(:, 5, Cnt) = arrayShrink(taskOpMotorMap(:, 1:size(allenMask,2)),allenMask,'merge');
            
            % isolate instructed movement events
            opMotorLabels = opMotorLabels(ismember(opMotorLabels,regLabels(unique(regIdx(~rejIdx))))); %reject variables that were not used in the full model
%             opLabels{iAnimals} = opMotorLabels;
            for iMoves = 1 : length(opMotorLabels)
                try
                    cIdx = find(opMotorIdx == find(ismember(opMotorLabels,opMotorLabels(iMoves))));
                    cData = opMotorR(:,cIdx(moveShift));
                    cData = reshape(cData,frames,[]);
                    opMoves{iAnimals}(:,Cnt,iMoves) = mean(cData,2);
                    opMoveTrialVar{iAnimals}(Cnt,iMoves,1) = var(sum(cData,1));
                    opMoveTrialVar{iAnimals}(Cnt,iMoves,2) = sum(cData(:))/size(cData,2);
                catch
                    opMoves{iAnimals}(:,Cnt,iMoves) = NaN(frames,1);
                    opMoveTrialVar{iAnimals}(Cnt,iMoves,:) = NaN;
                end
            end
            
            % isolate uninstructed movement events
            spontMotorLabels = spontMotorLabels(ismember(spontMotorLabels,regLabels(unique(regIdx(~rejIdx))))); %reject variables that were not used in the full model
%             spontLabels{iAnimals} = spontMotorLabels;
            for iMoves = 1 : length(spontMotorLabels)
                try
                    cIdx = find(spontMotorIdx == find(ismember(spontMotorLabels,spontMotorLabels(iMoves))));
                    cData = spontMotorR(:,cIdx(moveShift));
                    cData = reshape(cData,frames,[]);
                    spontMoves{iAnimals}(:,Cnt,iMoves,1) = mean(cData,2); %psth
                    spontMoveTrialVar{iAnimals}(Cnt,iMoves,1) = var(sum(cData,1)); %trial-by-trial variability
                    spontMoveTrialVar{iAnimals}(Cnt,iMoves,2) = sum(cData(:))/size(cData,2); %sum of all events, normalized by trialcount
                catch
                    spontMoves{iAnimals}(:,Cnt,iMoves) = NaN(frames,1);
                end
            end
            
            % get beta power - task model
            taskBeta = cat(3,taskBeta{:});
            taskBeta = nanmean(taskBeta,3);
            temp = nanmean(abs(U * taskBeta'),2);
            temp = alignAllenTransIm(arrayShrink(temp,mask,'split'), opts.transParams);
            taskPower{iAnimals}(:,Cnt) = arrayShrink(temp(:, 1:size(allenMask,2)),allenMask,'merge');
            
            % get beta power - instructed movement model
            opMotorBeta = cat(3,opMotorBeta{:});
            opMotorBeta = nanmean(opMotorBeta,3);
            temp = nanmean(abs(U * opMotorBeta'),2);
            temp = alignAllenTransIm(arrayShrink(temp,mask,'split'), opts.transParams);
            opMotorPower{iAnimals}(:,Cnt) = arrayShrink(temp(:, 1:size(allenMask,2)),allenMask,'merge');
            
            % get beta power - uninstructed movement model
            spontMotorBeta = cat(3,spontMotorBeta{:});
            spontMotorBeta = nanmean(spontMotorBeta,3);
            temp = nanmean(abs(U * spontMotorBeta'),2);
            temp = alignAllenTransIm(arrayShrink(temp,mask,'split'), opts.transParams);
            spontMotorPower{iAnimals}(:,Cnt) = arrayShrink(temp(:, 1:size(allenMask,2)),allenMask,'merge');
            
            % get variance in face camera
            a = NaN(1,length(faceM));
            for x = 1 : length(faceM)
                a(x) = nanmean(faceM{x});
            end
            allTwitch{iAnimals}(Cnt) = var(a);
            
            fprintf('Using recording %d/%d \n',iRecs,length(recs))
            toc
        end
    end
    recIdx{iAnimals} = recIdx{iAnimals}(useIdx);
end

%% look at results
titles = {'full model' 'movements' 'task-alone' 'uninstructed-alone' 'instructed-alone'};
lastSeg = 5;
lastDate = inf;
for x = 1:length(allDates) %get last recording to use based on shortest recording
    lastDate = min([lastDate min(allDates{x}(find(recIdx{x} == lastSeg, 1, 'last')))]);
    %     lastDate = min([lastDate min(allDates{x}(10))]);
end

%explained variance change
clear nMean xVals nSem
figure('renderer','painters');
aColors = {'g' 'b' 'k'};
runCnt = 0;
for y = [3 5 4]
    runCnt = runCnt +1;
    meanData = NaN(length(recIdx),100);
    for x = 1 : length(recIdx)
        ax1(runCnt) = subplot(1,4,runCnt); hold on; title(titles{y});
        %         cDates = allDates{x}(allDates{x} <= lastDate) - allDates{x}(1);
        cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
        plot(ax1(runCnt), cDates, allMaps{x}(allDates{x} <= lastDate,y),'linewidth',2,'color',[0.5 0.5 0.5]);
        axis square; ylim([0 0.5]);
        
        %         vline(cDates(find(recIdx{x} == 2, 1, 'first')),'--k');
        %         vline(cDates(find(recIdx{x} == 3, 1, 'first')),'--r');
%         vline(cDates(find(recIdx{x} == 4, 1, 'first')),'--g');
        %         vline(cDates(find(recIdx{x}== 5, 1, 'first')),'--c');
        xlabel('days');ylabel('cvR^2');
        meanData(x,cDates+1) = allMaps{x}(allDates{x} <= lastDate,y);
    end
    Cnt = 0; %compute average in time bins of 3 sessions
    for xx = 1 : 3 : size(meanData,2)-1
        Cnt = Cnt + 1;
        temp = meanData(:,xx:xx+2);
        nMean(Cnt) = nanmean(temp(:));
        nSem(Cnt) = sem(temp(:));
    end
    plot(ax1(runCnt),(find(~isnan(nMean))-1)*3,(nMean(~isnan(nMean))),'linewidth',4,'color','k');
    xlim(ax1(runCnt),[0 find(~isnan(nMean),1,'last')*3]); grid on;
    
    ax2(runCnt) = subplot(1,4,4); hold on; axis square
    cData = nanmean(meanData); cError = sem(meanData); xVals = find(~isnan(cData));
%     stdshade([],0.5,aColors{runCnt},xVals,5,cData(xVals),cError(xVals))
    stdshade([],0.5,aColors{runCnt},(find(~isnan(nMean))-1)*3,1,nMean(~isnan(nMean)),nSem(~isnan(nMean)))
end

% unique contributions task/movement
clear meanData cLines cTitles
figure('renderer','painters');
runCnt = 0;
for y = 1 : 7
    meanData{y} = NaN(length(recIdx),100);
    for x = 1 : length(recIdx)
        cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
        meanData{y}(x,cDates+1) = allMaps{x}(allDates{x} <= lastDate,y);
        if y == 3
            subplot(3,2,2); hold on
            cData = meanData{1}(x,:) - meanData{2}(x,:);
            plot(cDates, cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
            title(['Unique ' titles{y}]);
        elseif y == 6
            subplot(3,2,4); hold on
            cData = meanData{1}(x,:) - meanData{6}(x,:);
            plot(cDates, cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
            title(['Unique ' titles{5}]);
        elseif y == 7
            subplot(3,2,6); hold on
            cData = meanData{1}(x,:) - meanData{7}(x,:);
            plot(cDates, cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
            title(['Unique ' titles{4}]);
        end
    end
    
    if ismember(y,[1 3:5])
        runCnt = runCnt + 1;
        subplot(2,2,[1 3]);hold on
        Cnt = 0; %compute average in time bins of 3 sessions
        clear cData
        for xx = 1 : 3 : size(meanData{y},2)-1
            Cnt = Cnt + 1;
            temp = meanData{y}(:,xx:xx+2);
            cData(Cnt) = nanmean(temp(:));
        end
        cLines(runCnt) = plot((find(~isnan(cData))-1)*3,(cData(~isnan(cData))),'linewidth',2);
        cTitles{runCnt} = titles{y};
        ylim([0 0.6]); xlim([0 find(~isnan(cData),1,'last')*3]); title('Explained variance');
        ylabel('cvR^2');xlabel('days');
    end
    
    if y == 3
        subplot(3,2,2); hold on
        cData = nanmean(meanData{1} - meanData{2},1);
        plot(find(~isnan(cData))-1, smooth(cData(~isnan(cData))),'linewidth',4,'color','k'); axis square
        axis square;
        xlim([0 find(~isnan(cData),1,'last')]);
        ylabel('dvR^2'); xlabel('days');
    elseif y == 6
        subplot(3,2,4); hold on
        cData = nanmean(meanData{1} - meanData{6},1);
        plot(find(~isnan(cData))-1, smooth(cData(~isnan(cData))),'linewidth',4,'color','k');
        axis square; xlim([0 find(~isnan(cData),1,'last')]);
        ylabel('dvR^2');xlabel('days');
    elseif y == 7
        subplot(3,2,6); hold on
        cData = nanmean(meanData{1} - meanData{7},1);
        plot(find(~isnan(cData))-1, smooth(cData(~isnan(cData))),'linewidth',4,'color','k');
        axis square; xlim([0 find(~isnan(cData),1,'last')]);
        ylabel('dvR^2');xlabel('days');
    end
end
legend(cLines,cTitles); xlim([0 50]);


%% align against task performance
figure('renderer','painters');
for x = 1 : length(recIdx)
    for y = 1:2

        subplot(2,2,y); hold on;
        cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
        plot(cDates,smooth(allMaps{x}(allDates{x} <= lastDate,y+1)),'linewidth',2);
        axis square; xlim([-1 find(~isnan(meanData{1}(1,:)),1,'last')]);
        xlabel('days');ylabel('cvR^2'); ylim([0 0.6]);
%         if y == 1; title('unaligned task');
%         else; title('unaligned move'); ylim([0.55 0.75]); end

        subplot(2,2,y+2); hold on;
        cDates = allDates{x}(allDates{x} <= lastDate) - allDates{x}(find(recIdx{x}==4,1)); %shift to align with optimal performance
        plot(cDates,smooth(allMaps{x}(allDates{x} <= lastDate,y+1)),'linewidth',2);
        axis square;
        xlabel('days');ylabel('cvR^2'); ylim([0 0.6]);
%         if y == 1; title('aligned task'); ylim([0.3 0.55]);
%         else; title('aligned move'); ylim([0.55 0.75]); end
        vline(0);

    end
end


%% dimensionality and variance change
figure('renderer','painters');
meanData = NaN(2,length(recIdx),100);
for x = 1 : length(recIdx)
    subplot(2,1,1);
    hold on; title('widefield variance');
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
    plot(cDates, allVar{x}(allDates{x} <= lastDate),'linewidth',2,'color',[0.5 0.5 0.5]);
    axis square;
    xlabel('days');ylabel('cvR^2');
    meanData(1,x,cDates+1) = allVar{x}(allDates{x} <= lastDate);
    
    subplot(2,1,2);
    hold on; title('widefield dimensionality');
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
    plot(cDates, allDim{x}(allDates{x} <= lastDate),'linewidth',2,'color',[0.5 0.5 0.5]);
    axis square;
    xlabel('days');ylabel('cvR^2');
    meanData(2,x,cDates+1) = allDim{x}(allDates{x} <= lastDate);
end
subplot(2,1,1);
meanData = squeeze(nanmean(meanData,2));
plot(find(~isnan(meanData(1,:)))-1,smooth(meanData(1,~isnan(meanData(1,:)))),'linewidth',4,'color','k');
xlim([-1 find(~isnan(meanData(1,:)),1,'last')]);

subplot(2,1,2);
plot(find(~isnan(meanData(2,:)))-1,smooth(meanData(2,~isnan(meanData(2,:)))),'linewidth',4,'color','k');
xlim([-1 find(~isnan(meanData(2,:)),1,'last')]);


%% check for changes in trial-by-trial variability
figure('renderer','painters');
meanData = NaN(length(recIdx),100);
for x = 1 : length(recIdx)
    hold on; title('Trial-by-trial variance');
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
%     cData = mean(allData{x}(allDates{x} <= lastDate,:,2),2)';
    cData = allTrialVar{x}(allDates{x} <= lastDate);
    plot(cDates, cData,'linewidth',2,'color',[0.5 0.5 0.5]);
    axis square;
    xlabel('days');ylabel('Cross-trial variance');
    meanData(x,cDates+1) = cData;
    vline(cDates(find(recIdx{x} == 4, 1, 'first')),'--r');
end

meanData = squeeze(nanmean(meanData,1));
plot(find(~isnan(meanData(1,:)))-1,smooth(meanData(1,~isnan(meanData(1,:)))),'linewidth',4,'color','k');
xlim([-1 find(~isnan(meanData(1,:)),1,'last')]);


%% check for changes in task-overlap
figure('renderer','painters');
meanData = NaN(length(recIdx),100,2);
cMaps = NaN(sum(~allenMask(:)),100,2,length(recIdx),'single');
for x = 1 : length(recIdx)
    hold on; title('Task/uninstructed overlap');
    dInd = allDates{x} <= lastDate;
    cDates = allDates{x}(dInd) - firstDate(x);
    dInd(cDates <= 6) = false;
    
    cDates = allDates{x}(dInd) - firstDate(x);
    cData = allMaps{x}(dInd, 6) - allMaps{x}(dInd, 3); %spontMotor - task model. gives task-independent R^2.

    plot(cDates, -cData,'linewidth',2,'color','b');
    meanData(x,cDates+1,1) = -cData;

    cData = allMaps{x}(dInd, 4) - cData; %task-independent uninstructed minus all uninstructed. Gives task-aligned R^2.
    plot(cDates, cData,'linewidth',2,'color','c');
    meanData(x,cDates+1,2) = cData;

    axis square;
    xlabel('days');ylabel('Aligned/independent cvR^2');
    hline(0);
    
    % get maps
    cMaps(:,cDates,1,x) = squeeze(motorMaps{x}(:, 4, dInd)) - taskMaps{x}(:, dInd); %task-independent
    cMaps(:,cDates,2,x) = squeeze(motorMaps{x}(:, 4, dInd)) - cMaps(:,cDates,1,x); %task-aligned
end
meanData = squeeze(nanmean(meanData,1));
plot(find(~isnan(meanData(:,1)))-1,smooth(meanData(~isnan(meanData(:,1)),1)),'linewidth',4,'color','k');
plot(find(~isnan(meanData(:,1)))-1,smooth(meanData(~isnan(meanData(:,1)),2)),'linewidth',4,'color','k');

firstFrame = find(~isnan(nMean(250,250,:,1)),1);
lastFrame = find(~isnan(nMean(250,250,:,1)),1,'last');
cIdx = firstFrame-1 : floor((lastFrame - firstFrame) / 3) : lastFrame;
% show maps for change in task-indepdendent spontaneous movements
for y = 1 : 2
figure('renderer','painters');
cRange = [0 0.4]; 
subplot(1,3,1)
cData = nanmean(cMaps(:,cIdx(1)*3:cIdx(2)*3,y,:),4);
mapImg = imshow(arrayShrink(nanmean(cData,2),allenMask,'split'),cRange);
colormap(mapImg.Parent,inferno); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colorbar;

subplot(1,3,2)
cData = nanmean(cMaps(:,cIdx(end-1)*3:cIdx(end)*3,y,:),4);
mapImg = imshow(arrayShrink(nanmean(cData,2),allenMask,'split'),cRange);
colormap(mapImg.Parent,inferno); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colorbar;

subplot(1,3,3)
cData = nanmean(nanmean(cMaps(:,cIdx(end-1)*3:cIdx(end)*3,y,:),4),2) - nanmean(nanmean(cMaps(:,cIdx(1)*3:cIdx(2)*3,y,:),4),2);
mapImg = imshow(arrayShrink(cData,allenMask,'split'),[-0.2 0.2]);
colormap(mapImg.Parent,colormap_blueblackred); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colorbar;
end
    
%% combine predictive power maps into larger time bins
clear meanData nMean
meanData = NaN([size(allenMask),length(recIdx),100,3],'single');
firstFrame = inf; lastFrame = -inf;
for x = 1 : length(recIdx)
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
    %     firstFrame = min([firstFrame cDates(1)])+1;
    %     lastFrame = max([lastFrame cDates(end)]);
    cInd = allDates{x} <= lastDate;
    meanData(:,:,x,cDates+1,1) = arrayShrink(taskMaps{x}(:,cInd),allenMask,'split');
    meanData(:,:,x,cDates+1,2) = arrayShrink(motorMaps{x}(:,2,cInd),allenMask,'split');
    meanData(:,:,x,cDates+1,3) = arrayShrink(motorMaps{x}(:,3,cInd),allenMask,'split');
end
Cnt = 0; %compute average in time bins of 3 sessions
for xx = 1 : 3 : size(meanData,4)-1
    Cnt = Cnt + 1;
    temp = meanData(:,:,:,xx:xx+2,:);
    temp = reshape(temp,[size(temp,1),size(temp,2),size(temp,3)*size(temp,4),size(temp,5)]);
    nMean(:,:,Cnt,:) = nanmean(temp,3);
end
firstFrame = find(~isnan(nMean(250,250,:,1)),1);
lastFrame = find(~isnan(nMean(250,250,:,1)),1,'last');
cIdx = firstFrame-1 : floor((lastFrame - firstFrame) / 3) : lastFrame;
figure
titles = {'task cvR^2' 'uninstructed movement cvR^2' 'instructed movement cvR^2'};
Cnt = 0;
cRange = [0 0.4];
diffRange = [-0.2 0.2];
% for x = 1:size(nMean,4)
for x = [2 1 3] %this order is to get the mask from uninstructed movements first :)
    for y = [1 length(cIdx)-1]
        Cnt = Cnt + 1;
        subplot(size(nMean,4),3,Cnt)
        mapImg = imshow(nanmean(nMean(:,:, cIdx(y):cIdx(y+1),x),3),cRange);
        colormap(mapImg.Parent,inferno); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(titles{x}); colorbar;
        
        if y == length(cIdx)-1
            Cnt = Cnt + 1;
            subplot(size(nMean,4),3,Cnt); hold on;
            cMap = nanmean(nMean(:,:, cIdx(y):cIdx(y+1),x),3) - nanmean(nMean(:,:, cIdx(1):cIdx(1+1),x),3);
            mapImg = imshow(cMap,diffRange);
            colormap(mapImg.Parent,colormap_blueblackred); axis image
            set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
            title(['delta ' titles{x}]); colorbar;
            cMap = abs(cMap)>0.1;
            cMap = imerode(cMap,strel('disk',10));
            cMap = imdilate(cMap,strel('disk',10));
            areaInfo = regionprops(cMap, 'Area');
            areaInfo = [sort(cat(1,areaInfo.Area),'descend'); 2];
            cMap = bwareaopen(cMap,areaInfo(2)-1); %don't use more than 2 areas
            deltaAreas{x} = cMap;
            areaOutline = bwboundaries(deltaAreas{2});
            for xx = 1 : length(areaOutline)
                plot(areaOutline{xx}(:,2),areaOutline{xx}(:,1),'w')
            end
            deltaTraces{x,1} = arrayShrink(nMean(:,:,:,x), deltaAreas{2}, 'merge');
            deltaTraces{x,2} = arrayShrink(nMean(:,:,:,x), ~deltaAreas{2}, 'merge');
            deltaTraces{x,3} = arrayShrink(meanData(:,:,:,:,x), deltaAreas{2}, 'merge');
            deltaTraces{x,4} = arrayShrink(meanData(:,:,:,:,x), ~deltaAreas{2}, 'merge');
        end
    end
end

% plot area traces
figure
for y = 1:size(deltaTraces,1)
    
    ax1 = subplot(2,size(deltaTraces,1),y); hold on;
    ax2 = subplot(2,size(deltaTraces,1),y+size(deltaTraces,1)); hold on;
    for x = 1:size(meanData,3)
        cData = squeeze(nanmean(deltaTraces{y,3}(:,x,:),1));
        plot(ax1, find(~isnan(cData)),cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
        
        cData = squeeze(nanmean(deltaTraces{y,4}(:,x,:),1));
        plot(ax2, find(~isnan(cData)),cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
    end
    cData = squeeze(nanmean(deltaTraces{y,1},1));
    plot(ax1,((find(~isnan(cData))-1)*3)+1,(cData(~isnan(cData))),'linewidth',4,'color','k');
    ylim(ax1,[0 0.6]); axis(ax1,'square');
    title(ax1,[titles{y} ' - Area 1']);
    
    cData = squeeze(nanmean(deltaTraces{y,2},1));
    plot(ax2,((find(~isnan(cData))-1)*3)+1,(cData(~isnan(cData))),'linewidth',4,'color','k');
    ylim(ax2,[0 0.6]); axis(ax2,'square');
    title(ax2,[titles{y} ' - Area 2']);
end

% barplot to check for R^2 differences within and without selected areas
clear predMean predSem predBase predH predP predVals
for xx = 1:2
for x = 1 : 3
    Cnt = 0;
    for y = [1 length(cIdx)-1]
        cData = squeeze(nanmean(deltaTraces{x,xx+2}(:,:,cIdx(y)*3:cIdx(y+1)*3),1));
        cData = nanmean(cData,2);
        if y == 1
            predBase(:,xx,x) = cData;
        else
            Cnt = Cnt +1;
            predVals(:,x,xx,Cnt) = cData - predBase(:,xx,x);
            predMean(xx,x,Cnt) = mean(cData) - mean(predBase(:,xx,x));
            predSem(xx,x,Cnt) = sem(cData);
            [predH(xx,x,Cnt), predP(xx,x,Cnt)] = ttest(cData);
        end
    end
end
end

figure(99);
cIdx(1:2) = [1 6];
subplot(1,2,1); hold on;
errorbar((1:3)-0.15, predMean(1,:),predSem(1,:),'.k','linewidth',2);
errorbar((1:3)+0.15, predMean(2,:),predSem(2,:),'.k','linewidth',2);
bar(predMean');
for x = 1 : size(predVals,2)
    plot((1:3)-0.15, predVals(:,:,1) ,'o','color',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceColor','w');
    plot((1:3)+0.15, predVals(:,:,2) ,'o','color',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceColor','w');
end
axis square; hold off; title('Change in cvR^2 - Different areas');
legend({'Area1' 'Area2'});
set(gca,'xTick',1:3)
set(gca,'xTickLabel',titles)


%% combine unique contribution maps into larger time bins
clear meanData
meanData = NaN([size(allenMask),length(recIdx),100,3],'single');
firstFrame = inf; lastFrame = -inf;
for x = 1 : length(recIdx)
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
    firstFrame = min([firstFrame cDates(1)])+1;
    lastFrame = max([lastFrame cDates(end)]);
    cInd = allDates{x} <= lastDate;
    meanData(:,:,x,cDates+1,1) = arrayShrink(fullMaps{x}(:,cInd) - squeeze(motorMaps{x}(:,1,cInd)),allenMask,'split');
    meanData(:,:,x,cDates+1,2) = arrayShrink(fullMaps{x}(:,cInd) - squeeze(motorMaps{x}(:,5,cInd)),allenMask,'split');
    meanData(:,:,x,cDates+1,3) = arrayShrink(fullMaps{x}(:,cInd) - squeeze(motorMaps{x}(:,4,cInd)),allenMask,'split');
end
meanData = squeeze(nanmean(meanData,3));

figure
titles = {'task dR^2' 'uninstructed movement dR^2' 'instructed movement dR^2'};
Cnt = 0;
for x = 1:size(meanData,4)
    cIdx = firstFrame : floor((lastFrame - firstFrame) / 3) : lastFrame;
    if x == 1
        cRange = [0 0.02];
        diffRange = [-0.025 0.025];
    elseif x == 2
        cRange = [0 0.3];
        diffRange = [-0.2 0.2];
    elseif x == 3
        cRange = [0 0.075];
        diffRange = [-0.025 0.025];
    end
    for y = [1 3]
        Cnt = Cnt + 1;
        subplot(size(meanData,4),3,Cnt)
        mapImg = imshow(nanmean(meanData(:,:, cIdx(y):cIdx(y+1),x),3),cRange);
        colormap(mapImg.Parent,inferno); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(titles{x}); colorbar;
        
        if y == 3
            Cnt = Cnt + 1;
            subplot(size(meanData,4),3,Cnt)
            mapImg = imshow(nanmean(meanData(:,:, cIdx(y):cIdx(y+1),x),3) - nanmean(meanData(:,:, cIdx(1):cIdx(1+1),x),3),diffRange);
            colormap(mapImg.Parent,colormap_blueblackred); axis image
            set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
            title(['delta ' titles{x}]); colorbar;
        end
    end
end

%% look at changes in unique task contribution during learning
figure
cIdx = firstFrame : floor((lastFrame - firstFrame) / 4) : lastFrame;
cRange = [0 0.02];

for y = 1 : length(cIdx)-1
    
    subplot(1,length(cIdx)-1,y)
    mapImg = imshow(nanmean(meanData(:,:, cIdx(y):cIdx(y+1),1),3),cRange);
    colormap(mapImg.Parent,inferno); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    colorbar;
    
end

%% beta weight change
clear meanData nMean
meanData = NaN([size(allenMask),length(recIdx),100,3],'single');
firstFrame = inf; lastFrame = -inf;
for x = 1 : length(recIdx)
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
    %     firstFrame = min([firstFrame cDates(1)])+1;
    %     lastFrame = max([lastFrame cDates(end)]);
    cInd = allDates{x} <= lastDate;
    meanData(:,:,x,cDates+1,1) = arrayShrink(taskPower{x}(:,cInd),allenMask,'split');
    meanData(:,:,x,cDates+1,2) = arrayShrink(opMotorPower{x}(:,cInd),allenMask,'split');
    meanData(:,:,x,cDates+1,3) = arrayShrink(spontMotorPower{x}(:,cInd),allenMask,'split');
end
Cnt = 0; %compute average in time bins of 3 sessions
for xx = 1 : 3 : size(meanData,4)-1
    Cnt = Cnt + 1;
    temp = meanData(:,:,:,xx:xx+2,:);
    temp = reshape(temp,[size(temp,1),size(temp,2),size(temp,3)*size(temp,4),size(temp,5)]);
    nMean(:,:,Cnt,:) = nanmean(temp,3);
end

firstFrame = find(~isnan(nMean(250,250,:,1)),1);
lastFrame = find(~isnan(nMean(250,250,:,1)),1,'last');
cIdx = firstFrame-1 : floor((lastFrame - firstFrame) / 3) : lastFrame;
figure
titles = {'task cvR^2' 'uninstructed movement cvR^2' 'instructed movement cvR^2'};
Cnt = 0;
cRange = [0 0.0015];
diffRange = [-0.4 0.4];
for x = 1:size(nMean,4)
    for y = [1 length(cIdx)-1]
        Cnt = Cnt + 1;
        subplot(size(nMean,4),3,Cnt)
        mapImg = imshow(nanmean(nMean(:,:, cIdx(y):cIdx(y+1),x),3),cRange);
        colormap(mapImg.Parent,inferno); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(titles{x}); colorbar;
        
        if y == length(cIdx)-1
            Cnt = Cnt + 1;
            subplot(size(nMean,4),3,Cnt); hold on;
            firstMap = nanmean(nMean(:,:, cIdx(1):cIdx(1+1),x),3);
            cMap = nanmean(nMean(:,:, cIdx(y):cIdx(y+1),x),3);
            cMap = (cMap - firstMap)./firstMap;
            mapImg = imshow(cMap,diffRange);
            colormap(mapImg.Parent,colormap_blueblackred); axis image
            set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
            title(['delta ' titles{x}]); colorbar;
            
            contour(deltaAreas{2},'w');
            deltaTraces{x,1} = arrayShrink(nMean(:,:,:,x), deltaAreas{2}, 'merge');
            deltaTraces{x,2} = arrayShrink(nMean(:,:,:,x), ~deltaAreas{2}, 'merge');
            deltaTraces{x,3} = arrayShrink(meanData(:,:,:,:,x), deltaAreas{2}, 'merge');
            deltaTraces{x,4} = arrayShrink(meanData(:,:,:,:,x), ~deltaAreas{2}, 'merge');
        end
    end
end

%% plot area traces
figure; clear barData
for y = 1:size(deltaTraces,1)
    
    ax1 = subplot(2,size(deltaTraces,1),y); hold on;
    ax2 = subplot(2,size(deltaTraces,1),y+size(deltaTraces,1)); hold on;
    for x = 1:size(meanData,3)
        baseWeight = nanmean(nanmean(deltaTraces{y,3}(:,x,cIdx(1)*3:cIdx(2)*3),1));
        cData = (squeeze(nanmean(deltaTraces{y,3}(:,x,:),1))) ./ baseWeight; cData(cData > 2) = NaN;
        plot(ax1, find(~isnan(cData)),cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
        barData(x,y,1) = nanmean(cData(cIdx(end-1)*3:cIdx(end)*3))-1; %keep data for barplot
        
        baseWeight = nanmean(nanmean(deltaTraces{y,4}(:,x,cIdx(1)*3:cIdx(2)*3),1));
        cData = (squeeze(nanmean(deltaTraces{y,4}(:,x,:),1))) ./ baseWeight; cData(cData > 2) = NaN;
        plot(ax2, find(~isnan(cData)),cData(~isnan(cData)),'linewidth',2,'color',[0.5 0.5 0.5]);
        barData(x,y,2) = nanmean(cData(cIdx(end-1)*3:cIdx(end)*3))-1; %keep data for barplot
    end
    baseWeight = nanmean(nanmean(nanmean(deltaTraces{y,3}(:,:,cIdx(1)*3:cIdx(2)*3),1)));
    cData = squeeze(nanmean(nanmean(deltaTraces{y,3},1),2) ./ baseWeight);
    plot(ax1,((find(~isnan(cData))))+1,smooth(cData(~isnan(cData))),'linewidth',4,'color','k');
    ylim(ax1,[0 3]);
    axis(ax1,'square');
    title(ax1,[titles{y} ' - Area 1']);
    
    baseWeight = nanmean(nanmean(nanmean(deltaTraces{y,4}(:,:,cIdx(1)*3:cIdx(2)*3),1)));
    cData = squeeze(nanmean(nanmean(deltaTraces{y,4},1),2) ./ baseWeight);
    plot(ax2,((find(~isnan(cData))))+1,smooth(cData(~isnan(cData))),'linewidth',4,'color','k');
    ylim(ax2,[0 3]);
    axis(ax2,'square');
    title(ax2,[titles{y} ' - Area 2']);
end

%% barplot to check for beta weight differences within and without selected areas
clear weightMean weightSem baseWeight weightH weightP weightVals

for xx = 1:2
    for x = 1 : 3
        baseWeight = squeeze(nanmean(nanmean(deltaTraces{x,xx+2}(:,:,cIdx(1)*3:cIdx(2)*3),1),3))';
        cData = squeeze(nanmean(deltaTraces{x,xx+2}(:,:,cIdx(end-1)*3:cIdx(end)*3),1));
        cData = bsxfun(@rdivide,cData, baseWeight); cData(cData > 2) = NaN;
        cData = nanmean(cData,2)-1;
        
        weightVals(:,x,xx) = cData;
        weightMean(xx,x) = mean(cData);
        weightSem(xx,x) = sem(cData);
        [weightH(xx,x), weightP(xx,x,Cnt)] = ttest(cData);
    end
end

figure(99);
subplot(1,2,2); cla; hold on;
errorbar((1:3)-0.15, weightMean(1,:),weightSem(1,:),'.k','linewidth',2);
errorbar((1:3)+0.15, weightMean(2,:),weightSem(2,:),'.k','linewidth',2);
bar(weightMean');
for x = 1 : size(weightVals,2)
    plot((1:3)-0.15, weightVals(:,:,1) ,'o','color',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceColor','w');
    plot((1:3)+0.15, weightVals(:,:,2) ,'o','color',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceColor','w');
end
axis square; hold off; title('Change in cvR^2 - Different areas');
legend({'Area1' 'Area2'});
set(gca,'xTick',1:3)
set(gca,'xTickLabel',titles)

