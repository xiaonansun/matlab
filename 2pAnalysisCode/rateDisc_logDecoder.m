function [cvAcc, bMaps, mdlAll, trialCnt, cvAccShuf] = rateDisc_logDecoder(dataIn, U, bhv, useTrials, targMod, regType, stepSize, decType, learnType, shufReps)
% 2021-06-23 added trialNumbers
% run logistic regression decoder on neural data.

% cvAcc and bMaps have 5 decoder outputs each: 1 = all trials choice, 2 =
% balanced trials choice, 3 = balanced trials stimulus, 4 = correct trials
% choice, 5 = error trials choice
% 2021-07-08 Richard changed input variable "data" to "dataIn" to prevent
% variable name conflicts during debugging

% INPUTS
% dataIn: input neural data. 3 dimensional matrix. neuron_ID x time x
% trial_ID (rows are neurons/cells, columns are samples in time, 3rd
% dimension is trials)

% U: empty for now

% bhv: behavior data, in this case Bpod struct

% useTrials: use only selected trials. In this case, trials will be
% selected by the function rateDisc_equalizeTrials, which balances trials.

% targMod: target modality: 0 for all trials, number >= 1 for other
% modalities.

% regType: regularization type, lasso vs ridge regression

% stepSize: subsampling of recordings with high sampling rates. for 2p data
% (usually recorded at 30 fps), the stepSize input variable can be left
% empty, subsequently defaulting the value to 1, where data at every time
% point is used for analysis. If stepSize > 1, then time points will be
% skipped based on the stepSize

% decType: decoder type. Specify a task variable for the decoder to work on.
% The options are 'allChoice', 'allStim', 'correct', 'error', 'choice', 'preChoice',
% 'nextChoice'.

% learnType: This variable is the input for the fitclinear 'Learner' option
% options are 'logistic' versus 'svm'.

% shufReps: number of repetitions for shuffling 

if ~exist('decType','var') || isempty(decType)
    decType = 'allChoice';
end

if ~exist('stepSize','var') || isempty(stepSize)
    stepSize = 1;
end

if ~exist('learnType','var') || isempty(learnType)
    learnType = 'logistic';
end

if ~exist('shufReps','var') || isempty(shufReps)
    shufReps = 50;
end

dSize = floor(size(dataIn,2) / stepSize); %maximum number of bins

if isempty(U)
    bMaps = NaN(size(dataIn,1),dSize, 'single');
    bMapsShuf = NaN(size(dataIn,1),dSize, 'single');
%     allAUC = NaN(size(dataIn,1),dSize, 'single');
else
    bMaps = NaN(size(U,1), dSize, 'single'); %beta maps for different decoders
    bMapsShuf = NaN(size(U,1), dSize, 'single'); %beta maps for different decoders
end
cvAcc = NaN(1, dSize, 'single'); cvAccShuf = NaN(shufReps, dSize, 'single');
trialCnt = NaN(1, dSize, 'single'); 
% betaNeuron = NaN(1, dSize, 'single'); % Added by Richard 2021-07-06, for single-neuron beta weights
mdlAll = cell(1,dSize); % Added by Richard 2021-07-06, for single-neuron models
mdlAllShuf = cell(shufReps,dSize); % Added by Richard 2021-07-06, for single-neuron models
% leftIdxShuf = nan(shufReps,dSize); 

%% get some indices, normalize data and check trialcount for each time point
if targMod == 0
    expIdx = true(1, length(bhv.StimType)); %use all trials
else
    expIdx = bhv.StimType == targMod; %trials with target modality
end
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

for i = 1:shufReps
    leftIdxShuf(i,:) = leftIdx(randperm(length(leftIdx)));
end

%% run decoder
rng(1) % for reproducibility
Cnt = 0;
for iSteps = 1 : stepSize : dSize*stepSize
%%    
    if size(dataIn,1)==1
        cData = mean(dataIn(:, iSteps : iSteps + stepSize - 1, :),2,'omitnan');
    else
        cData = squeeze(mean(dataIn(:, iSteps : iSteps + stepSize - 1, :),2,'omitnan'));
    end
    
    Cnt = Cnt + 1;

    if strcmpi(decType, 'allChoice')
        %% all trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, [], useTrials); %equalize L/R choices
        
%         indices{iSteps} = choiceIdx; %keeping track of which trials were selected in the previous line
        trialCnt(iSteps) = sum(choiceIdx); %nr of trials in choice decoder
        if sum(choiceIdx) >= useTrials
            
            X = cData(:, choiceIdx);
            Xmean = mean(X,2);
            X = bsxfun(@minus,X,Xmean);
            Xstd = std(X,0,2);
            useIdx = Xstd>0;
            X = X(useIdx,:);
            Xstd = Xstd(useIdx);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(leftIdx(choiceIdx)');

            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
            cvAcc(Cnt) = 1-kfoldLoss(Mdl);
            mdlAll{Cnt} = Mdl; % added by Richard 2021-07-06
            
            if shufReps > 0
                parfor j = 1:shufReps
                    choiceIdxShuf = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdxShuf(j,:), [], useTrials); %equalize L/R choices
                    if sum(choiceIdxShuf) <= sum(choiceIdx)
                        XShuf = cData(:, choiceIdxShuf);
                        XShufmean = mean(XShuf,2);
                        XShuf = bsxfun(@minus,XShuf,XShufmean);
                        XShufstd = std(XShuf,0,2);
                        useIdx = XShufstd>0;
                        XShuf = XShuf(useIdx,:);
                        XShufstd = XShufstd(useIdx);
                        XShuf = bsxfun(@rdivide,XShuf,XShufstd);
                    else
                        XShuf = X;
                    end
                    YShuf = double(leftIdxShuf(choiceIdxShuf)'); % Shuffled
                    MdlShuf = fitclinear(XShuf, YShuf, 'ObservationsIn','columns', 'kfold', 10, 'Regularization', regType, 'Learner', learnType);
                    cvAccShuf(j,Cnt) = 1-kfoldLoss(MdlShuf);
                    %                 disp(['Shuffling trials, iteration #' num2str(j)]);
                    %                 mdlAllShuf{Cnt} = MdlShuf; % added by Richard 2021-07-06
                end
            end
            if ~isempty(U)
                clear a
                for x = 1:length(Mdl.Trained) % cycles through individual cross-validation trainings
                    a(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd); % divides beta by the standard deviation
                end
                bMaps(:, Cnt) = nanmean(a,2);
            else
                clear a aShuf
                for x = 1:length(Mdl.Trained)
                    a(:,x) = bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
                end
                bMaps(useIdx, Cnt) = nanmean(a,2); % computes the mean beta across all (in this case, 10) cross-validation runs
                
%                 if shufReps > 0
%                     for x = 1:length(MdlShuf.Trained)
%                         aShuf(:,x) = bsxfun(@rdivide,MdlShuf.Trained{x}.Beta, XShufstd);
%                     end
%                     bMapsShuf(useIdx, Cnt) = nanmean(aShuf,2); % computes the mean beta across all (in this case, 10) cross-validation runs
%                 end
            end
        end
        
        %%
    elseif strcmpi(decType, 'allStim')
        % all trials - stim decoder
        stimIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), targIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(iSteps) = sum(stimIdx); %nr of balanced trials in stim decoder
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

%%
    elseif strcmpi(decType, 'correct')
        % correct trials - choice decoder
        choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)) & corrIdx, leftIdx, [], useTrials); %equalize L/R choices and correct/incorrect trials
        trialCnt(iSteps) = sum(choiceIdx); %nr of balanced trials in choice decoder
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
        trialCnt(iSteps) = sum(choiceIdx); %nr of balanced trials in choice decoder
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
        trialCnt(iSteps) = sum(choiceIdx); %nr of balanced trials in choice decoder
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
        trialCnt(iSteps) = sum(stimIdx); %nr of balanced trials in stim decoder
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
        trialCnt(iSteps) = sum(pChoiceIdx); %nr of balanced trials in choice decoder
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
        trialCnt(iSteps) = sum(nChoiceIdx); %nr of balanced trials in choice decoder
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
