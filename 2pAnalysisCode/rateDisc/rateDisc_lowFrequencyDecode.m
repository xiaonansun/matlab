
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
baseFrames = 1:60; %baseline frames - use first 2 seconds
freqRange = 3:6; %frequency range for power analysis
decSteps = 3; %width of frequency band for decoder in Hz
maxFreq = 10; %maximal frequency for freq ranges

%% do things
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
currSteps = 0 : decSteps : maxFreq; %steps for decoder frequency bands

% get animal info
dataOverview = rateDiscRecordings;
animals = dataOverview(:,1);
animals = animals(1:4);

decAcc = cell(1, length(animals));
trialCnt = cell(1, length(animals));
absMaps = cell(1, length(animals));
betaMaps = cell(1, length(animals));

for iAnimals = 1 : length(animals)
    
    % get recs
    load([cPath animals{iAnimals} filesep 'blockData' filesep 'trialInfo.mat'],'recs')
    absMaps{iAnimals} = NaN(sum(~allenMask(:)), length(recs), length(currSteps) - 1, 2, 'single');
    betaMaps{iAnimals} = NaN(sum(~allenMask(:)), length(recs), length(currSteps) - 1, 'single');
    decAcc{iAnimals} = NaN(length(recs), length(currSteps) - 1);
    trialCnt{iAnimals} = NaN(1,length(recs));
    
    for iRecs = 1 : length(recs)
        
        fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; disp(fPath); %session data path
        cFile = dir([fPath '*_SpatialDisc*.mat']);
        load(fullfile(fPath,cFile.name));
        load([fPath 'Vc.mat'],'Vc','U','bTrials');
        load([fPath 'opts2.mat'],'opts');
        
        U = alignAllenTransIm(U,opts.transParams);
        U = arrayShrink(U(1:size(allenMask,1),1:size(allenMask,2),:),allenMask); %merge pixels
        
        Vc = Vc(:, baseFrames,:);
        bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dat
        
        Cnt = zeros(1,2);
        allFFT = cell(1,size(Vc,3));
        useIdx = true(1, size(Vc,3));
        for iTrials = 1 : size(Vc,3)
            cIdx = ~isnan(squeeze(nanmean(Vc(:,:,iTrials)))); %non-NaN frames
            
            if sum(cIdx) >= 30
                % do FFT
                allFFT{iTrials} = fft(squeeze(Vc(:,cIdx,iTrials)),[],2);
                
                % check if trial was rewarded
                rIdx = bhv.Rewarded(iTrials) + 1;
                Cnt(rIdx) = Cnt(rIdx) + 1;
                
                currSteps = 0 : decSteps : maxFreq; %steps for decoder frequency bands
                currSteps = round(currSteps * (sum(cIdx) / 30)) + 1;
                currSteps(1) = 2; %don't use first bin
                
                % collect maps for different frequencies
                for iSteps = 1 : length(currSteps) - 1
                    freqIdx = currSteps(iSteps) : currSteps(iSteps+1);
                    temp = U * allFFT{iTrials}(:, freqIdx); %get correct part from fft
                    temp = sum(abs(temp),2); %make single amplitude map for selected frequency range
                    absMaps{iAnimals}(:, iRecs, iSteps, rIdx) = runMean(absMaps{iAnimals}(:, iRecs, iSteps, rIdx), temp, Cnt(rIdx)); %running average
                end
            else
                useIdx(iTrials) = false; %don't use this trial
            end
        end
        
        % check if any frequency range in the baseline is predictive of success
        corrIdx = rateDisc_equalizeTrials(useIdx, bhv.Rewarded); %equalize success/error trials
        trialCnt{iAnimals}(iRecs) = sum(corrIdx); %nr of trials in current decoder
        
        for iSteps = 1 : length(currSteps) - 1
            % move through frequency bands
            Cnt = 0;
            X = zeros(size(Vc,1), sum(corrIdx));
            for iTrials = find(corrIdx)
                
                cIdx = ~isnan(squeeze(nanmean(Vc(:,:,iTrials)))); %non-NaN frames
                if sum(cIdx) >= 30
                    currSteps = 0 : decSteps : maxFreq; %steps for decoder frequency bands
                    currSteps = round(currSteps * (sum(cIdx) / 30)) + 1;
                    currSteps(1) = 2; %don't use first bin
                    
                    Cnt = Cnt + 1;
                    X(:,Cnt) = sum(abs(allFFT{iTrials}(:,currSteps(iSteps):currSteps(iSteps+1))),2);
                end
            end
            Xmean = mean(X,2);
            X = bsxfun(@minus,X,Xmean);
            Xstd = std(X,0,2);
            X = bsxfun(@rdivide,X,Xstd);
            Y = double(bhv.Rewarded(corrIdx)');
            
            % get cross-validated model to compute predictive power
            Mdl = fitclinear(X, Y, 'ObservationsIn','columns', 'kfold', 5, 'Learner', 'svm', 'Regularization','ridge');
            decAcc{iAnimals}(iRecs, iSteps) = 1 - kfoldLoss(Mdl);
            clear temp
            for x = 1:length(Mdl.Trained)
                temp(:,x) = U * bsxfun(@rdivide,Mdl.Trained{x}.Beta, Xstd);
            end
            betaMaps{iAnimals}(:, iRecs, iSteps) = nanmean(temp,2);
        end
    end
end
