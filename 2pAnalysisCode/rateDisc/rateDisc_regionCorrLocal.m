function [regCorr, dimCorr, B] = rateDisc_regionCorrLocal(A, C, areas, regs, dimCnt)
%% compute cross-region correlation
%interleave dims from both hemisspheres
cRegs =  [regs; -regs + 2^8];
cRegs = cRegs(:)';

%make index to only use 'dimCnt' highest dims in each area
targIdx= [];
sourceIdx= [];
Cnt = 0;
for x = cRegs
    y = find(areas == x); %combine dims from both hemisspheres
    targIdx = [targIdx Cnt + (1:min([length(y) dimCnt]))]; %index for target array  that has 'dimCnt' rows for each region and HS
    sourceIdx = [sourceIdx y(1:min([length(y) dimCnt]))]; %index for dimensions in A to match regions
    Cnt = Cnt + dimCnt;
end

%% get correlation matrices and sort into larger output matrix
B = corrcoef(C(sourceIdx,:)'); %correlation between dimensions
dimCorr = NaN(length(cRegs) * dimCnt, 'single'); %sorted correlation map between dimensions
dimCorr(targIdx,targIdx) = B;

%% predict different regions
regCorr = NaN(length(regs), 'single'); %correlation between regions
for iRegs = 1 : length(regs)
    
    cIdx = areas == iRegs | -(areas-2^8) == iRegs;
    cData = C(cIdx,:)'; %current dimensions. Use those to predict other regions.
    cData = bsxfun(@minus,cData,nanmean(cData,1)); %subtract mean
    
    for oRegs = 1 : length(regs) %predict other regions
        cIdx = areas == oRegs | -(areas-2^8) == oRegs;
        targData = C(cIdx,:)'; %current dimensions. Use those to predict other regions.
        targData = bsxfun(@minus,targData,nanmean(targData,1)); %subtract mean
        
        [~, dimBeta] = ridgeMML(targData, cData, true); %get ridge penalties and beta weights.
        Vm = cData * dimBeta; %modeled data
        [~, varP1, ~, covP] = rateDisc_modelCorr(targData', Vm', A(:,:,cIdx)); %compute variance and covariance
        regCorr(iRegs, oRegs) = sum(real(covP)) / sum(varP1); %total amount of explained variance in percent (similar to R^2)
    end
end
end
