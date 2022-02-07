function [allAUC, shufAUC, sAUC] = twoP_aucAnalysisNew(dataIn, bhv, shuffleAUC, doSDU, animal, session)
% 2021-06-23 added trialNumbers
% run logistic regression decoder on widefield data.
% cvAcc and bMaps have 5 decoder outputs each: 1 = all trials choice, 2 =
% balanced trials choice, 3 = balanced trials stimulus, 4 = correct trials
% choice, 5 = error trials choice
% 2021-07-08 Richard changed input variable "data" to "dataIn" to prevent
% variable name conflicts during debugging

S = twoP_settings;
imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;

if ~exist('decType','var') || isempty(decType)
    decType = 'allChoice';
end

if ~exist('stepSize','var') || isempty(stepSize)
    stepSize = 1;
end

if ~exist('learnType','var') || isempty(learnType)
    learnType = 'logistic';
end

if ~exist('shuffleAUC','var') || isempty(shuffleAUC)
    shuffleAUC = 0;
elseif shuffleAUC == 1
    nShuf = 100;
elseif shuffleAUC ==0
    shufAUC = [];
end


dSize = floor(size(dataIn,2) / stepSize); %maximum number of bins


    bMaps = NaN(size(dataIn,1),dSize, 'single');
    allAUC = NaN(size(dataIn,1),dSize, 'single');

cvAcc = NaN(1, dSize, 'single');
trialCnt = NaN(1, dSize, 'single');
betaNeuron = NaN(1, dSize, 'single'); % Added by Richard 2021-07-06, for single-neuron beta weights
mdlAll = cell(1,dSize); % Added by Richard 2021-07-06, for single-neuron models

%% get some indices, normalize data and check trialcount for each time point
    expIdx = true(1, length(bhv.StimType)); %use all trials

perfIdx = ~bhv.DidNotChoose & ~bhv.DidNotLever & expIdx; %performed trials with target modality
% perfIdx = intersect(perfIdx,trialNumbers); % added by Richard 2021-06-23
dataIn = dataIn(:, :, perfIdx); %reject non-performed trials

leftIdx = (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left
pLeftIdx = [false leftIdx(1:end-1)];
pLeftIdx = pLeftIdx(perfIdx); %previous choices
nLeftIdx = [leftIdx(2:end) false];
nLeftIdx = nLeftIdx(perfIdx); %future choices
leftIdx = leftIdx(perfIdx); %current choices
targIdx = bhv.CorrectSide(perfIdx) == 1; %trials were left was correct side
corrIdx = bhv.Rewarded(perfIdx); %rewarded trials

if shuffleAUC == 1
    for iShuf = 1:nShuf
        leftIdxShuf(:,iShuf) = leftIdx(randperm(length(leftIdx)))';
    end
end
%% run decoder
rng(1) % for reproducibility
% iSteps = 0;

%%
for iSteps = 1 : stepSize : dSize*stepSize
% for iSteps = 50
    if size(dataIn,1)==1
        cData = nanmean(dataIn(:, iSteps : iSteps + stepSize - 1, :),2);
    else
        cData = squeeze(nanmean(dataIn(:, iSteps : iSteps + stepSize - 1, :),2));
    end
    
%     iSteps = iSteps + 1;

    if strcmpi(decType, 'allChoice')
        % all trials - choice decoder
%         choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, [], useTrials); %equalize L/R choices
%         trialiSteps(iSteps) = sum(choiceIdx); %nr of trials in choice decoder
%         if sum(choiceIdx) >= useTrials
            
%             X = cData(:, choiceIdx);
            X = cData(:, :);
            
            % compute SDU
            if doSDU == 1
                Xmean = mean(X,2,'omitnan');
                X = bsxfun(@minus,X,Xmean);
                Xstd = std(X,0,2,'omitnan');
                useIdx = Xstd>0;
                X = X(useIdx,:);
                Xstd = Xstd(useIdx);
                X = bsxfun(@rdivide,X,Xstd);
            elseif doSDU == 0 
%                 Xmean = mean(X,2,'omitnan');
%                 X = bsxfun(@minus,X,Xmean);
                Xstd = std(X,0,2,'omitnan');
                useIdx = Xstd>0;
                X = X(useIdx,:);
%                 Xstd = Xstd(useIdx);
%                 X = bsxfun(@rdivide,X,Xstd);
            end
            
            idxAUC = ~isnan(mean(X));
            Xauc = X(:,idxAUC);
%             Y = double(leftIdx(choiceIdx)');
            Y = double(leftIdx');
            Yauc = Y(idxAUC);
            sAUC(iSteps).X = Xauc; sAUC(iSteps).Y = Yauc;
            % auc analysis
%             allAUC(useIdx, iSteps) = colAUC(X',Y);
            allAUC(useIdx, iSteps) = colAUC(Xauc',Yauc);
            if shuffleAUC == 1
                Yshuf = leftIdxShuf(idxAUC,:);
                for iShuf = 1:nShuf
                    shufAUC(useIdx, iSteps, iShuf) = colAUC(Xauc',Yshuf(:,iShuf));
                end
            end
            
    elseif strcmpi(decType, 'allStim')
        % all trials - stim decoder
        stimIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), targIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(stimIdx); %nr of balanced trials in stim decoder
        if sum(stimIdx) >= useTrials
            
            X = cData(:, stimIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(targIdx(stimIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'correct')
        % correct trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)) & corrIdx, leftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(choiceIdx); %nr of balanced trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'error')
        % error trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)) & ~corrIdx, leftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(choiceIdx); %nr of balanced trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'choice')
        % balanced trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, corrIdx, useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(choiceIdx); %nr of balanced trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'stim')
        % balanced trials - stim decoder
        stimIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), targIdx, corrIdx, useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(stimIdx); %nr of balanced trials in stim decoder
        if sum(stimIdx) >= useTrials
            
            X = cData(:, stimIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(targIdx(stimIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'preChoice')
        % balanced trials - previous choice decoder
        pChoiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), pLeftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(pChoiceIdx); %nr of balanced trials in choice decoder
        if sum(pChoiceIdx) >= useTrials
            
            X = cData(:, pChoiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(pLeftIdx(pChoiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'nextChoice')
        % balanced trials - previous choice decoder
        nChoiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), nLeftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialiSteps(iSteps) = sum(nChoiceIdx); %nr of balanced trials in choice decoder
        if sum(nChoiceIdx) >= useTrials
            
            X = cData(:, nChoiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(nLeftIdx(nChoiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(iSteps) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, iSteps) = nanmean(a,2);
            end
        end
    end
end


if (~exist('animal','var') || isempty(animal)) && (~exist('session','var') || isempty(session))
    disp('AUC data will not be saved due to missing inputs: animal and session IDs.');
else
        AUC.allAUC=allAUC; AUC.shufAUC=shufAUC; AUC.sAUC = sAUC;
        S = twoP_settings;
        imagingRootDir = S.dir.imagingRootDir;
        imagingSubDir = S.dir.imagingSubDir;
        saveDir= fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir);
        save(fullfile(saveDir,['AUC_' num2str(size(AUC.allAUC,2)) '.mat']),'AUC');
        disp(['AUC data saved as: ' fullfile(saveDir,['AUC_' num2str(size(AUC.allAUC,2)) '.mat'])]);
end

end


