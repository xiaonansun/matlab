function [cvAcc, bMaps, trialCnt] = rateDisc_logDecoder(data, U, bhv, useTrials, targMod, regType, stepSize, decType)
% run logistic regression decoder on widefield data.
% cvAcc and bMaps have 5 decoder outputs each: 1 = all trials choice, 2 =
% balanced trials choice, 3 = balanced trials stimulus, 4 = correct trials
% choice, 5 = error trials choice

if ~exist('decType','var') || isempty(decType)
    decType = 'allChoice';
end

if ~exist('stepSize','var') || isempty(stepSize)
    stepSize = 1;
end
dSize = floor(size(data,2) / stepSize); %maximum number of bins

if isempty(U)
    bMaps = NaN;
else
    bMaps = NaN(size(U,1), dSize, 'single'); %beta maps for different decoders
end
cvAcc = NaN(1, dSize, 'single');
trialCnt = NaN(1, dSize, 'single');
learnType = 'logistic';

%% get some indices, normalize data and check trialcount for each time point
if targMod == 0
    expIdx = true(1, length(bhv.StimType)); %use all trials
else
    expIdx = bhv.StimType == targMod; %trials with target modality
end
perfIdx = ~bhv.DidNotChoose & ~bhv.DidNotLever & expIdx; %performed trials with target modality
data = data(:, :, perfIdx); %reject non-performed trials

leftIdx = (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left
pLeftIdx = [false leftIdx(1:end-1)];
pLeftIdx = pLeftIdx(perfIdx); %previous choices
nLeftIdx = [leftIdx(2:end) false];
nLeftIdx = nLeftIdx(perfIdx); %future choices
leftIdx = leftIdx(perfIdx); %current choices
targIdx = bhv.CorrectSide(perfIdx) == 1; %trials were left was correct side
corrIdx = bhv.Rewarded(perfIdx); %rewarded trials

%% run decoder
rng(1) % for reproducibility
Cnt = 0;
for iSteps = 1 : stepSize : dSize*stepSize

    cData = squeeze(nanmean(data(:, iSteps : iSteps + stepSize - 1, :),2));
    Cnt = Cnt + 1;

    if strcmpi(decType, 'allChoice')
        % all trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, [], useTrials); %equalize L/R choices
        trialCnt(Cnt) = sum(choiceIdx); %nr of trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            Xmean = mean(X,2);
            X = bsxfun(@minus,X,Xmean);
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                clear a
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'allStim')
        % all trials - stim decoder
        stimIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), targIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(stimIdx); %nr of balanced trials in stim decoder
        if sum(stimIdx) >= useTrials
            
            X = cData(:, stimIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(targIdx(stimIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'correct')
        % correct trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)) & corrIdx, leftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(choiceIdx); %nr of balanced trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'error')
        % error trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)) & ~corrIdx, leftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(choiceIdx); %nr of balanced trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'choice')
        % balanced trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, corrIdx, useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(choiceIdx); %nr of balanced trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'stim')
        % balanced trials - stim decoder
        stimIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), targIdx, corrIdx, useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(stimIdx); %nr of balanced trials in stim decoder
        if sum(stimIdx) >= useTrials
            
            X = cData(:, stimIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(targIdx(stimIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'preChoice')
        % balanced trials - previous choice decoder
        pChoiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), pLeftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(pChoiceIdx); %nr of balanced trials in choice decoder
        if sum(pChoiceIdx) >= useTrials
            
            X = cData(:, pChoiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(pLeftIdx(pChoiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
        
    elseif strcmpi(decType, 'nextChoice')
        % balanced trials - previous choice decoder
        nChoiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), nLeftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(Cnt) = sum(nChoiceIdx); %nr of balanced trials in choice decoder
        if sum(nChoiceIdx) >= useTrials
            
            X = cData(:, nChoiceIdx);
            X = bsxfun(@minus,X,mean(X,2));
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(nLeftIdx(nChoiceIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            if ~isempty(U)
                for x = 1:length(Mdl.Trained)
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(:, Cnt) = nanmean(a,2);
            end
        end
    end
end
end