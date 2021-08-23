[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, ~, fPath, trainDates] = rateDiscRecordings;
cPath = 'Y:\data\BpodImager\Animals\';
stepSize = 15;
regType = 'ridge'; %lasso or ridge
animals = dataOverview(:,1);
recs = dataOverview(:,3);

% Cnt = 0; brokenRec = [];
for iAnimals = 1 : length(animals)
    %% get recordings and sort by date
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep];
    recs = ls(fPath);
    recs = recs(~ismember(recs(:,1), '.'), :);
    cDate = datenum(recs(:, 1:11)); %get dates from recordings
   
    if size(recs,2) > 11 %of there are multiple recordings from the same date
        for iRecs = find(ismember(recs(:,12), '_'))
           cDate(iRecs) = cDate(iRecs) + (str2double(recs(iRecs,13))/100);
        end
    end
    [cDate,ind] = sort(cDate,'ascend'); % sort by oldest first
    recs = recs(ind,:); %adjust order of filenames to get it to be chronological
    recs = mat2cell(recs, ones(1,size(recs,1)), size(recs,2)); %convert to cells
    
    if iAnimals == 1
        % get allen maps on first run
        load('allenDorsalMapSM.mat')
        mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
        [x1, y1] = size(mask);
        load([fPath recs{1} filesep 'snapshot_1.mat'])
        load([fPath recs{1} filesep 'opts2.mat'])
        snap = alignAllenTransIm(single(snap),opts.transParams);
        [x2, y2] = size(snap);
        mask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
        segIdxRealign = {1:54 55:79 80:98 113:131 132:162 163:179}; %segment index after realignment on both handle and stimulus
        rightHs = find(ismember(dorsalMaps.sidesSplit,'L')); %index for labels on the left HS
    end
        
    for iRecs = 1:length(recs)
        %% predict choice using re-aligned widefield
        cFile = dir([fPath recs{iRecs} filesep animals{iAnimals} '_SpatialDisc*.mat']);
        load([fPath recs{iRecs} filesep strtrim(cFile.name)]); %load behavior data
        load([fPath recs{iRecs} filesep 'Vc.mat'], 'bTrials', 'U', 'Vc');
        load([fPath recs{iRecs} filesep 'regData.mat']);
        load([fPath recs{iRecs} filesep 'opts2.mat']);
        
        U = alignAllenTransIm(U,opts.transParams);
        U = arrayShrink(U(1:size(mask,1),1:size(mask,2), :),mask);
        
        [alignIdx, trialIdx, frames] = Widefield_getRealignment(fullR, idx, recIdx, trialIdx, recLabels, frames);
        Vc = Vc(:,alignIdx);
        Vc = bsxfun(@minus,Vc,mean(Vc,2));
        
        Vc = reshape(Vc, size(Vc,1), frames, []);
        [cpRaw{iAnimals}, betaRaw{iAnimals}] = delayDec_logRegress(SessionData, Vc, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
    
%     %load motor model
%     load([fPath 'motorBeta.mat'], 'motorBeta');
%     load([fPath 'motorregData.mat'], 'motorR');
%     motorBeta = mean(cat(3,motorBeta{:}),3);
%     Vm = (motorR * motorBeta)';
%     Vm = Vm(:,alignIdx);
%     Vm = bsxfun(@minus,Vm,mean(Vm,2));
%     Vm = reshape(Vm, size(Vm,1), frames, []);
%     Vm = smoothCol(Vm,[],[],[],2);
%     Vm = reshape(Vm, size(Vm,1), []);
% 
    %load spontaneous motor model
    load([fPath 'orgVspontMotor.mat'], 'Vm');
    Vm = Vm(:,alignIdx);
    Vm = bsxfun(@minus,Vm,mean(Vm,2));
    Vm = reshape(Vm, size(Vm,1), frames, []);
    Vm = smoothCol(Vm,[],[],[],2);
    Vm = reshape(Vm, size(Vm,1), []);
    
    % get regression coefficients
    cData = smoothCol(reshape(Vc, size(Vc,1), frames, []),[],[],[],2);
    cData = reshape(cData, size(cData,1), []);
    regSp = nansum(Vm.*cData, 2) ./ nansum(Vm.*Vm, 2);
    spontCorrect = reshape(cData - Vm .* regSp, size(cData,1), frames, []);
    Vm = reshape(Vm, size(Vm,1), frames, []);

    [cpSpontCorrectTask{iAnimals}, betaSpontCorrectTask{iAnimals}] = delayDec_logRegress(SessionData, spontCorrect, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
    [cpSpontMotor{iAnimals}, betaSpontMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
    clear cData Vm Vspont motorCorrect spontCorrect

%     %subtract motor prediction from raw data
%     [cpTask{iAnimals}, betaTask{iAnimals}] = delayDec_logRegress(SessionData, motorCorrect, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
%     %subtract spontaneous motor prediction from raw data
%     [cpSpontCorrectTask{iAnimals}, betaSpontCorrectTask{iAnimals}] = delayDec_logRegress(SessionData, spontCorrect, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
%     clear cData Vm Vspont motorCorrect spontCorrect
    
%     %load non-choice model and test task reconstruction
%     load([fPath 'noChoicedimBeta.mat'], 'dimBeta');
%     load([fPath 'noChoiceregData.mat'], 'fullR', 'recLabels', 'idx', 'recIdx');
%     
%     cInd = ismember(recIdx(~idx), find(~ismember(recLabels,motorLabels(~ismember(motorLabels,opLabels))))); %find task + op movements
%     rVtask = (fullR(:, cInd) * dimBeta(cInd,:))';
%     rVtask = rVtask(:,alignIdx);
%     rVtask = reshape(rVtask, size(rVtask,1), frames, []);
%     [cpRecTask{iAnimals}, betaRecTask{iAnimals}] = delayDec_logRegress(SessionData, rVtask, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
%    
%     cInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels(~ismember(motorLabels,opLabels))))); %find task + op movements
%     rVspont = (fullR(:, cInd) * dimBeta(cInd,:))';
%     rVspont = rVspont(:,alignIdx);
%     rVspont = reshape(rVspont, size(rVspont,1), frames, []);
%     [cpRecSpont{iAnimals}, betaRecSpont{iAnimals}] = delayDec_logRegress(SessionData, rVspont, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
%    
    
%     %load motor model
%     load([fPath 'motorBeta.mat'], 'motorBeta');
%     load([fPath 'motorregData.mat'], 'motorR');
%     motorBeta = mean(cat(3,motorBeta{:}),3);
%     Vm = (motorR * motorBeta)';
%     Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
%     
%     [cpMotor{iAnimals}, betaMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
%     [cpTask{iAnimals}, betaTask{iAnimals}] = delayDec_logRegress(SessionData, Vc-Vm, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
    
%      %load operant motor model
%     load([fPath 'opMotorBeta.mat'], 'opMotorBeta');
%     load([fPath 'opMotorregData.mat'], 'opMotorR');
%     motorBeta = mean(cat(3,opMotorBeta{:}),3);
%     Vm = (opMotorR * motorBeta)';
%     Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
%     [cpOpMotor{iAnimals}, betaOpMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);
%      
%     %load spont motor model
%     load([fPath 'spontMotorBeta.mat'], 'spontMotorBeta');
%     load([fPath 'spontMotorregData.mat'], 'spontMotorR');
%     motorBeta = mean(cat(3,spontMotorBeta{:}),3);
%     Vm = (spontMotorR * motorBeta)';
%     Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
%     [cpSpontMotor{iAnimals}, betaSpontMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);

    %     %load no choice/no spont movement task model
%     load([fPath 'TaskNoChoiceNoSpontdimBeta.mat'], 'dimBeta');
%     load([fPath 'TaskNoChoiceNoSpontregData.mat'], 'fullR');
%     Vm = (fullR * dimBeta)';
%     Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
%     [cpNoSpontTask{iAnimals}, betaNoSpontTask{iAnimals}] = delayDec_logRegress(SessionData, Vm, U, bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, regType);

    
    %load video only model
    load([fPath filesep 'orgregData.mat'], 'fullR', 'recLabels', 'recIdx', 'idx');
    cIdx = ismember(recIdx(~idx), find(ismember(recLabels,{'Move' 'bhvVideo'}))); %get index for task regressors
    videoR = fullR(:, cIdx)';
    videoR = reshape(videoR(:,alignIdx), size(videoR,1), frames, []);
    cpVideo{iAnimals} = delayDec_logRegress(SessionData, videoR, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize, 'ridge', 'svm');

    toc
end


%% predict choice figure
figure
cvChoice1 = cat(3,cpRaw{:});
subplot(2,2,1)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Raw flourescence'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

% cvChoice1 = cat(3,cpRecSpont{:});
% subplot(2,2,2)
% lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
% lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
% lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
% ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
% vline([55 81 99 114 132 162].*(1/30));
% title('10x cvChoice - Spont. corrected raw'); legend(lines,{'All' 'Expert' 'Novice'})
% xlabel('Time(s)'); ylabel('Classifier accuracy');
% clear lines

cvChoice1 = cat(3,cpSpontMotor{:});
subplot(2,2,3)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Movement corrected raw'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

cvChoice1 = cat(3,cpSpontCorrectTask{:});
subplot(2,2,4)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Inst. movement reconstruction'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

%% predict novice with expert trials and vice-versa
cvChoice1 = cat(3,cpRaw{:});
figure
lines(1) = stdshade(squeeze(cvChoice1(:,4,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,5,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30))
title('Predicted choice - Cross-modal prediction'); legend(lines,{'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');

%% show results for video prediction
figure
cvChoice1 = cat(3,cpVideo{:});
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Video prediction'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

%% look at choice maps
clear temp1
temp = nanmean(cat(4,betaRaw{:}),4);
temp = smoothCol(temp,3,'gauss',[],2);
temp = arrayShrink(temp,mask,'split');
% temp1 = temp;
temp1(:,1:586/2,:,:) = (temp(:,1:586/2,:,:) - temp(:,end : -1 : 586/2 + 1,:,:));
% temp1(:,end : -1 : 586/2 + 1,:,:) = (temp(:,end : -1 : 586/2 + 1,:,:) - temp(:,1:586/2,:,:));

temp2 = abs(temp1(:,:,:,2)) - abs(temp1(:,:,:,3)); %difference between expert and novices
figure;
for iMod = 1:3
    subplot(1,3,iMod);
    mapImg = imagesc(smooth2a(nanmean(temp2(:, :, round((segIdxRealign{iMod + 2})/stepSize)),3), 5, 5)); axis image;
    hold on; caxis([0 0.02]);
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData));
end
colormap inferno

% compareMovie(temp1(:,:,:,1:3),'outline',dorsalMaps.edgeOutlineSplit(ismember(dorsalMaps.sidesSplit,'L')))