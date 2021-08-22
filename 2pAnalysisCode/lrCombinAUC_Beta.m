function a = lrCombinAUC_Beta(a,LR)

%% Combine AUC and Beta

%AUC

a.allAUC.analysisType = 'AUC';

a.allAUC.deep.NR = []; a.allAUC.deep.R = []; a.allAUC.deep.name = 'Deep';
try
for i = 1:length(a.iDeep)
    idx1 = LR(a.iDeep(i)).iEpoch(1,3); idx2 = LR(a.iDeep(i)).iEpoch(2,3);
    a.allAUC.deep.NR = [a.allAUC.deep.NR; mean(LR(a.iDeep(i)).allAUC(LR(a.iDeep(i)).idx_notredcell,idx1:idx2),2)];
    a.allAUC.deep.R = [a.allAUC.deep.R; mean(LR(a.iDeep(i)).allAUC(LR(a.iDeep(i)).idx_redcell,idx1:idx2),2)];
end
end

a.allAUC.int.NR = []; a.allAUC.int.R = []; a.allAUC.int.name = 'Intermediate';
try
for i = 1:length(a.iInt)
    idx1 = LR(a.iInt(i)).iEpoch(1,3); idx2 = LR(a.iInt(i)).iEpoch(2,3);
    a.allAUC.int.NR = [a.allAUC.int.NR; mean(LR(a.iInt(i)).allAUC(LR(a.iInt(i)).idx_notredcell,idx1:idx2),2)];
    a.allAUC.int.R = [a.allAUC.int.R; mean(LR(a.iInt(i)).allAUC(LR(a.iInt(i)).idx_redcell,idx1:idx2),2)];
end
end

a.allAUC.sup.NR = []; a.allAUC.sup.R = []; a.allAUC.sup.name = 'Superficial';
try
for i = 1:length(a.iSup)
    idx1 = LR(a.iSup(i)).iEpoch(1,3); idx2 = LR(a.iSup(i)).iEpoch(2,3);
    a.allAUC.sup.NR = [a.allAUC.sup.NR; mean(LR(a.iSup(i)).allAUC(LR(a.iSup(i)).idx_notredcell,idx1:idx2),2)];
    a.allAUC.sup.R = [a.allAUC.sup.R; mean(LR(a.iSup(i)).allAUC(LR(a.iSup(i)).idx_redcell,idx1:idx2),2)];
end
end
% Beta
a.dBeta.analysisType = 'Beta';

a.dBeta.deep.NR = []; a.dBeta.deep.R = []; a.dBeta.deep.name = 'Deep';
try
for i = 1:length(a.iDeep)
    idx1 = LR(a.iDeep(i)).iEpoch(1,3); idx2 = LR(a.iDeep(i)).iEpoch(2,3);
    a.dBeta.deep.NR = [a.dBeta.deep.NR; mean(LR(a.iDeep(i)).bMaps(LR(a.iDeep(i)).idx_notredcell,idx1:idx2),2)];
    a.dBeta.deep.R = [a.dBeta.deep.R; mean(LR(a.iDeep(i)).bMaps(LR(a.iDeep(i)).idx_redcell,idx1:idx2),2)];
end
end

a.dBeta.int.NR = []; a.dBeta.int.R = []; a.dBeta.int.name = 'Intermediate';
try
for i = 1:length(a.iInt)
    idx1 = LR(a.iInt(i)).iEpoch(1,3); idx2 = LR(a.iInt(i)).iEpoch(2,3);
    a.dBeta.int.NR = [a.dBeta.int.NR; mean(LR(a.iInt(i)).bMaps(LR(a.iInt(i)).idx_notredcell,idx1:idx2),2)];
    a.dBeta.int.R = [a.dBeta.int.R; mean(LR(a.iInt(i)).bMaps(LR(a.iInt(i)).idx_redcell,idx1:idx2),2)];
end
end

a.dBeta.sup.NR = []; a.dBeta.sup.R = []; a.dBeta.sup.name = 'Superficial';
try
for i = 1:length(a.iSup)
    idx1 = LR(a.iSup(i)).iEpoch(1,3); idx2 = LR(a.iSup(i)).iEpoch(2,3);
    a.dBeta.sup.NR = [a.dBeta.sup.NR; mean(LR(a.iSup(i)).bMaps(LR(a.iSup(i)).idx_notredcell,idx1:idx2),2)];
    a.dBeta.sup.R = [a.dBeta.sup.R; mean(LR(a.iSup(i)).bMaps(LR(a.iSup(i)).idx_redcell,idx1:idx2),2)];
end
end