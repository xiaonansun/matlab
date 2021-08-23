function [regCorr, dimCorr, B, regIdx] = rateDisc_regionCorr(A, C, areas, reg, allDims, dimCnt, useAreas, dorsalMaps)

if ~exist('dorsalMaps','var') || isempty(dorsalMaps)
    load('allenDorsalMapSM.mat');
end

%% compute cross-region correlation
%make index to only use 'dimCnt' highest dims in each area
targIdx = false(1, length(allDims) * dimCnt); %index for correlation values in dimcorr matrix
areaIdx =[];
Cnt = 0;
for x = allDims'
    y = find(areas == x);
    areaIdx = [areaIdx y(1:min([length(y) dimCnt]))];
    targIdx(Cnt + (1:min([length(y) dimCnt]))) = true;
    Cnt = Cnt + dimCnt;
end

%% get correlation matrices - sorted with 10 dims for each area
B = corrcoef(C(areaIdx,:)'); %correlation between dimensions
dimCorr = NaN(length(allDims) * dimCnt, 'single'); %sorted correlation map between dimensions
Cnt = 0;
for iRegs = allDims'
    y = areas(areaIdx) == iRegs;
    dimCorr(Cnt + (1:min([sum(y) dimCnt])), targIdx) = B(y,:);
    Cnt = Cnt + dimCnt;
end

%%
regCorr = NaN(length(reg), 'single'); %correlation between regions
for iRegs = 1 : length(reg)
    
    cIdx = reg{iRegs}'; %current areas
    cIdx = cIdx(ismember(cIdx, useAreas)); %only use areas within field of view
    regIdx{iRegs} = find(ismember(areas(areaIdx), cIdx)); %keep region index for correlation matrix
    cIdx = ismember(areas, cIdx); %get dimensions in C
    
    cData = C(cIdx,:)'; %current dimensions. Use those to predict other regions.
    cData = bsxfun(@minus,cData,nanmean(cData,1)); %subtract mean
    for oRegs = 1 : length(reg) %predict other regions
        
        cIdx = reg{oRegs}'; %current areas
        cIdx = cIdx(ismember(cIdx, useAreas)); %only use areas within field of view
        cIdx = ismember(areas, cIdx); %get dimensions in C

        targData = C(cIdx,:)'; %current dimensions. Use those to predict other regions.
        targData = bsxfun(@minus,targData,nanmean(targData,1)); %subtract mean
   
        [~, dimBeta] = ridgeMML(targData, cData, true); %get ridge penalties and beta weights.
        Vm = cData * dimBeta; %modeled data
        [~, varP1, ~, covP] = rateDisc_modelCorr(targData', Vm', A(:,:,cIdx)); %compute variance and covariance
        regCorr(iRegs, oRegs) = sum(covP) / sum(varP1); %total amount of explained variance in percent (similar to R^2)
        
    end
end   
end
