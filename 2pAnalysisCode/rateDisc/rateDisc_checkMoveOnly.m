% function rateDisc_checkMoveOnly
%code to asses how much neural variance is explained by movements variables in the linear encoding model across ssessions

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
moveShift = 16; %timepoint 0 in binary movement design matrx
spontMoveLabels = {'piezo' 'whisk' 'nose' 'body' 'face'};

% this is to isolate beta weights
[~, motorLabels] = delayDecRecordings;
allOp = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'}; %all operant motor regressors
allSpont = motorLabels(~ismember(motorLabels,allOp)); %all spontaneous motor regressors

%% get list of recordings and reference image and load data from all sessions
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
animals = {'mSM63' 'mSM64' 'mSM65' 'mSM66'};
% animals = {'Fez7' 'Fez10'};
% animals = {'Plex01' 'Plex02'};

opMoves = cell(1, length(animals));
opMoveTrialVar = cell(1, length(animals));
spontMoves = cell(1, length(animals));
spontMoveTrialVar = cell(1, length(animals));
spontAnaMoves = cell(1, length(animals));
spontAnaMoveTrialVar = cell(1, length(animals));

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
    for iRecs = 1 : length(recs)
        tic
        % get model data
        fPath = [aPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        gfPath = [gPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        checker = false;
        
        try
            
            load([fPath 'BehaviorVideo' filesep 'FilteredPupil.mat'],'bodyM');
            load([fPath 'interpVc.mat'],'frames');
            load([fPath 'Vc.mat'],'bTrials');            
            load([fPath 'regData.mat'],'rejIdx','regIdx','regLabels','trialIdx');
            load([fPath 'spontMotorregData.mat'],'spontMotorR','spontMotorIdx','spontMotorLabels');
            load([fPath 'opMotorregData.mat'],'opMotorR','opMotorIdx','opMotorLabels');
                        
            bhvFile = dir([fPath filesep animal '_SpatialDisc*.mat']);
            bhvFile = strtrim(bhvFile.name);
            load([fPath bhvFile]); %load behavior data
            trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
            bhv = selectBehaviorTrials(SessionData,bTrials(trialIdx)); %only use completed trials that are in the Vc dataset
            leftIdx = bhv.Rewarded & (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left and got rewarded

            if (size(opMotorR,1) / frames) < 50 || frames ~= 150; error; end
            checker = true;
        end
        if checker
            Cnt = Cnt +1; useIdx(iRecs) = true;
        
           % isolate instructed movement events
            cLabels = opMotorLabels(ismember(opMotorLabels,regLabels(unique(regIdx(~rejIdx))))); %reject variables that were not used in the full model
            for iMoves = 1 : length(opMoveLabels)
                try
                    cIdx = find(opMotorIdx == find(ismember(cLabels,opMoveLabels(iMoves))));
                    cData = opMotorR(:,cIdx(moveShift));
                    cData = reshape(cData,frames,[]);
                    cData = cData(:, leftIdx);
                    opMoves{iAnimals}(:,Cnt,iMoves,1) = mean(cData,2); %mean for time points
                    opMoves{iAnimals}(:,Cnt,iMoves,2) = var(cData,[],2); %variance in time points
                    opMoveTrialVar{iAnimals}(Cnt,iMoves,1) = var(sum(cData,1)); %variance across trials
                    opMoveTrialVar{iAnimals}(Cnt,iMoves,2) = mean(sum(cData,1)); %avg nr of events per trial
                catch
                    opMoves{iAnimals}(:,Cnt,iMoves,:) = NaN(frames,2);
                    opMoveTrialVar{iAnimals}(Cnt,iMoves,:) = NaN(1,2);
                end
            end
            
            % isolate uninstructed movement events
            cLabels = spontMotorLabels(ismember(spontMotorLabels,regLabels(unique(regIdx(~rejIdx))))); %reject variables that were not used in the full model
            for iMoves = 1 : length(spontMoveLabels)
                try
                    cIdx = find(spontMotorIdx == find(ismember(cLabels,spontMoveLabels(iMoves))));
                    cLength = length(cIdx);
                    if strcmpi(spontMoveLabels(iMoves),{'piezo'})
                        cRange = [1 moveShift];
                    elseif strcmpi(spontMoveLabels(iMoves),{'whisk'})
                        cRange = [1 ceil((cLength-1)/2) + moveShift]; %analog and high amp events
                    elseif strcmpi(spontMoveLabels(iMoves),{'nose'})
                        cRange = [1 ceil((cLength-1)/2) + moveShift]; %analog and high amp events
                    elseif strcmpi(spontMoveLabels(iMoves),{'body'})
                        cRange = [1 ceil((cLength-1)/2) + moveShift]; %analog and high amp events
                    elseif strcmpi(spontMoveLabels(iMoves),{'face'})
                        cRange = [1 ceil((cLength-1)/2) + moveShift]; %analog and high amp events
                    end
                    cData = spontMotorR(:,cIdx(cRange));
                    cData = reshape(cData,frames,[], size(cData,2));
                    cData = cData(:, leftIdx, :);

                    % get analog regressors
                    spontAnaMoves{iAnimals}(:,Cnt,iMoves,1) = mean(cData(:,:,1),2); %mean for time points
                    spontAnaMoves{iAnimals}(:,Cnt,iMoves,2) = var(cData(:,:,1),[],2); %variance in time points
                    spontAnaMoveTrialVar{iAnimals}(Cnt,iMoves,1) = var(sum(cData(:,:,1),1)); %variance across trials
                    spontAnaMoveTrialVar{iAnimals}(Cnt,iMoves,2) = mean(mean(cData(:,:,1),1)); %avg activity per trial
                    
                    % get binary regressors
                    spontMoves{iAnimals}(:,Cnt,iMoves,1) = mean(cData(:,:,2),2); %mean for time points
                    spontMoves{iAnimals}(:,Cnt,iMoves,2) = var(cData(:,:,2),[],2); %variance in time points
                    spontMoveTrialVar{iAnimals}(Cnt,iMoves,1) = var(sum(cData(:,:,2),1)); %variance across trials
                    spontMoveTrialVar{iAnimals}(Cnt,iMoves,2) = mean(sum(cData(:,:,2),1)); %avg nr of events per trial
                    
                catch
                    spontMoves{iAnimals}(:,Cnt,iMoves,:) = NaN(frames,2);
                    spontAnaMoves{iAnimals}(:,Cnt,iMoves,:) = NaN(frames,2);
                    spontMoveTrialVar{iAnimals}(Cnt,iMoves,:) = NaN(1,2);
                    spontAnaMoveTrialVar{iAnimals}(Cnt,iMoves,:) = NaN(1,2);
                end
            end
            allDates{iAnimals}(Cnt) = datenum(recs(iRecs).name);

            
            % get variance in face camera
            a = NaN(1,length(bodyM));
            for x = 1 : length(bodyM)
                a(x) = nanmean(bodyM{x});
            end
            allTwitch{iAnimals}(Cnt) = var(a);
            
            fprintf('Using recording %d/%d \n',iRecs,length(recs))
            toc
        end
    end
    recIdx{iAnimals} = recIdx{iAnimals}(useIdx);
end

%% look at results
lastSeg = 5;
lastDate = inf;
for x = 1:length(allDates) %get last recording to use based on shortest recording
    lastDate = min([lastDate min(allDates{x}(find(recIdx{x} == lastSeg, 1, 'last')))]);
end

%% movement change - face velocity
figure('renderer','painters');
meanData = NaN(length(recIdx),100);
for x = 1 : length(recIdx)
    hold on; title('movement change - Face camera');
    cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
    plot(cDates, allTwitch{x}(allDates{x} <= lastDate),'linewidth',2,'color',[0.5 0.5 0.5]);
    axis square;
    xlabel('days');ylabel('Standard deviation');
    meanData(x,cDates+1) = allTwitch{x}(allDates{x} <= lastDate);
end
meanData = nanmean(meanData,1);
plot(find(~isnan(meanData))-1,smooth(meanData(~isnan(meanData))),'linewidth',4,'color','k');
xlim([-1 find(~isnan(meanData),1,'last')]);

%% trial-by-trial variance - spontaneous movements
clear meanData
figure(50); set(gcf,'renderer','painters')
figure(51); set(gcf,'renderer','painters')
for y = 1 : length(spontMoveLabels)
    meanData = NaN(100,length(recIdx));
    for x = 1:length(recIdx)
        
        figure(50);
        subplot(ceil(length(spontMoveLabels) / 2),2,y); hold on;
        cData = spontMoveTrialVar{x}(:,y,1);
%         cData = mean(spontAnaMoves{x}(:,:,y,1),1)';
        cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
        meanData(cDates,x) = cData(allDates{x} <= lastDate); %trial-by-trial variance
        
        plot(cDates, meanData(cDates,x),'linewidth',2,'color',[0.5 0.5 0.5]);
        axis square; xlabel('days');ylabel('Mean event rate');
    end
    meanData = nanmean(meanData,2);
    plot(find(~isnan(meanData))-1,smooth(meanData(~isnan(meanData))),'linewidth',4,'color','k');
    xlim([-1 find(~isnan(meanData),1,'last')]); hold off
    title(spontMoveLabels{y})
    
    figure(51);
    plot(find(~isnan(meanData))-1,smooth(meanData(~isnan(meanData))),'linewidth',4); hold on
    axis square;
    
%     Cnt = 0; %compute average in time bins of 3 sessions
%     for xx = 1 : 3 : size(meanData,1)-1
%         Cnt = Cnt + 1;
%         temp = meanData(xx:xx+2,:);
%         nMean(Cnt) = nanmean(temp(:));
%     end
%     plot((find(~isnan(nMean))-1)*3,nMean(~isnan(nMean)),'linewidth',4,'color','k');
%     xlim([-1 find(~isnan(nMean),1,'last')*3]); 
%     title(spontMoveLabels{y})
end

%% trial-by-trial variance - operant movements
clear meanData
figure(55); set(gcf,'renderer','painters')
figure(56); set(gcf,'renderer','painters')
for y = 1 : length(opMoveLabels)
    meanData = NaN(100,length(recIdx));
    for x = 1:length(recIdx)
        
        figure(55);
        subplot(ceil(length(opMoveLabels) / 2),2,y); hold on;
        cData = opMoveTrialVar{x}(:,y,1);
        cDates = allDates{x}(allDates{x} <= lastDate) - firstDate(x);
        meanData(cDates,x) = cData(allDates{x} <= lastDate); %trial-by-trial variance
        
        plot(cDates, meanData(cDates,x),'linewidth',2,'color',[0.5 0.5 0.5]);
        axis square; xlabel('days');ylabel('Mean event rate');
    end
    meanData = nanmean(meanData,2);
    plot(find(~isnan(meanData))-1,smooth(meanData(~isnan(meanData))),'linewidth',4,'color','k');
    xlim([-1 find(~isnan(meanData),1,'last')]); 
    title(opMoveLabels{y})
    
    figure(56);
    plot(find(~isnan(meanData))-1,smooth(meanData(~isnan(meanData))),'linewidth',4); hold on
    axis square;
end

%% example spontaneous movement - show PSTH in each timepoint for different sessions
cVar = 'whisk';
iAnimal = 1;
cData = spontMoves{iAnimal}(:,:,ismember(spontMoveLabels,cVar),1);
cData = cData(3:end-2,:,:); %slight crop at the edges
cDates = allDates{iAnimal}(allDates{iAnimal} <= lastDate) - firstDate(iAnimal);
cIdx = cDates(1)-1 : floor((cDates(end) - cDates(1)) / 4) : cDates(end);

figure('renderer','painters');
ax = subplot(2,2,1);
rateDisc_linePlot(ax,cData, cDates, cIdx(1:2));
axis square; title([cVar ': Early sessions']);
ylabel('Mean event rate'); ylim([0 0.35]);
xlabel('Frames'); xlim([0 size(cData,1)]);
vline(59);

ax = subplot(2,2,2);
rateDisc_linePlot(ax,cData, cDates, cIdx(3:4))
axis square;  title([cVar ': Late sessions']);
ylabel('Mean event rate'); ylim([0 0.35]);
xlabel('Frames'); xlim([0 size(cData,1)]);
vline(59);

% example instructed movement - show variance in each timepoint for different sessions
cVar = 'lLick';
cData = opMoves{iAnimal}(:,:,ismember(opMoveLabels,cVar),1);
cData = cData(3:end-2,:,:); %slight crop at the edges
cDates = allDates{iAnimal}(allDates{iAnimal} <= lastDate) - firstDate(iAnimal);
cIdx = cDates(1)-1 : floor((cDates(end) - cDates(1)) / 4) : cDates(end);

ax = subplot(2,2,3);
rateDisc_linePlot(ax,cData, cDates, cIdx(1:2));
axis square; title([cVar ': Early sessions']);
ylabel('Mean event rate'); ylim([0 0.35]);
xlabel('Frames'); xlim([0 size(cData,1)]);
vline(59);

ax = subplot(2,2,4);
rateDisc_linePlot(ax,cData, cDates, cIdx(3:4))
axis square; title([cVar ': Late sessions']);
ylabel('Mean event rate'); ylim([0 0.35]);
xlabel('Frames'); xlim([0 size(cData,1)]);
vline(59);

%% example spontaneous movement - show fano factor in each timepoint for different sessions
cVar = 'whisk';
iAnimal = 1;
cData = spontMoves{iAnimal}(:,:,ismember(spontMoveLabels,cVar),2) ./ spontMoves{iAnimal}(:,:,ismember(spontMoveLabels,cVar),1);
cData = cData(3:end-2,:,:); %slight crop at the edges
cDates = allDates{iAnimal}(allDates{iAnimal} <= lastDate) - firstDate(iAnimal);
cIdx = cDates(1)-1 : floor((cDates(end) - cDates(1)) / 4) : cDates(end);

figure('renderer','painters');
ax = subplot(2,2,1);
rateDisc_linePlot(ax,cData, cDates, cIdx(1:2));
axis square; title([cVar ': Early sessions']);
ylabel('Fano factor'); ylim([0.6 1.01]);
xlabel('Frames'); xlim([0 size(cData,1)]);

ax = subplot(2,2,2);
rateDisc_linePlot(ax,cData, cDates, cIdx(3:4))
axis square;  title([cVar ': Late sessions']);
ylabel('Fano factor'); ylim([0.6 1.01]);
xlabel('Frames'); xlim([0 size(cData,1)]);

% example instructed movement - show variance in each timepoint for different sessions
cVar = 'lLick';
cData = opMoves{iAnimal}(:,:,ismember(opMoveLabels,cVar),2) ./ opMoves{iAnimal}(:,:,ismember(opMoveLabels,cVar),1);
cData = cData(3:end-2,:,:); %slight crop at the edges
cDates = allDates{iAnimal}(allDates{iAnimal} <= lastDate) - firstDate(iAnimal);
cIdx = cDates(1)-1 : floor((cDates(end) - cDates(1)) / 4) : cDates(end);

ax = subplot(2,2,3);
rateDisc_linePlot(ax,cData, cDates, cIdx(1:2));
axis square; title([cVar ': Early sessions']);
ylabel('Fano factor'); ylim([0.6 1.01]);
xlabel('Frames'); xlim([0 size(cData,1)]);

ax = subplot(2,2,4);
rateDisc_linePlot(ax,cData, cDates, cIdx(3:4))
axis square; title([cVar ': Late sessions']);
ylabel('Fano factor'); ylim([0.6 1.01]);
xlabel('Frames'); xlim([0 size(cData,1)]);
