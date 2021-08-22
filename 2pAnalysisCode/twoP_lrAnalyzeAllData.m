%% Analysis of combined logistic regression results

clear LR A AUC
close all;

if ~exist('LR','var')
    LR=twoP_lrLoadAllSessions(1);
end

animal = {'Fez51';'Fez57';'Fez59'};
% Combine AUC and Beta
% [A,AUC] = twoP_lrCombinAUC_Beta(LR,animal); 

% saveDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
% save(fullfile(saveDir,'combinedData',[datestr(now,'yyyy-mm-dd_HHMMSS') '_logisticRegression.mat']),'LR');


% plot AUC distributions
[A,AUC] = twoP_lrCombinAUC_Beta(LR,animal);

twoP_plotAUCHistogram(A,AUC.MM);

% AUC: fraction choice selective and magnitude of means

twoP_plotAUCFraction(A,AUC.MM)

%% AUC - plot CDF

[h,a.AllAUC.pSI,ks2stat]=kstest2(a.allAUC.sup.NR,a.allAUC.int.NR);
[h,a.AllAUC.pSD,ks2stat]=kstest2(a.allAUC.sup.NR,a.allAUC.deep.NR);
[h,a.AllAUC.pID,ks2stat]=kstest2(a.allAUC.int.NR,a.allAUC.deep.NR);

hAUCDepth = figure(1); hold on; 
hAUC(1) = cdfplot(a.allAUC.sup.NR);
hAUC(2) = cdfplot(a.allAUC.int.NR); 
hAUC(3) = cdfplot(a.allAUC.deep.NR);
set(hAUC,'Linewidth',2)
set(gca,'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02],...
  'xlim',[0 1]);
grid off;
txtSup=['Superficial (< 200 \mum)' newline 'n=' num2str(length(a.allAUC.sup.NR))]; 
txtInt = ['Intermediate (250-300 \mum)' newline 'n=' num2str(length(a.allAUC.int.NR))]; 
txtDeep = ['Deep (>500 \mum)' newline 'n=' num2str(length(a.allAUC.deep.NR))];
lAUC = legend(txtSup,...
    txtInt,...
    txtDeep,...
    'Box','off',...
    'Color','none',...
    'Location','Southeast');

% annotation('textbox',[0.2 0.5 0.3 0.3],'String',['Sup vs Int: ' num2str(pSI) 'Sup vs Deep: ' num2str(pSD) 'Int vs Deep: ' num2str(pID)]);
strAUC = ['Sup vs Int: p=' num2str(round(a.AllAUC.pSI,2,'significant')),...
    '\nSup vs Deep: p=' num2str(round(a.AllAUC.pSD,2,'significant')),...
    '\nInt vs Deep: p=' num2str(round(a.AllAUC.pID,2,'significant'))];
text(0.01,0.9,sprintf(strAUC));

title([a.allAUC.analysisType ' of non-' a.Setting.animal ' neurons: comparison across depths'])
xlabel(a.allAUC.analysisType);

exportgraphics(gcf,fullfile(a.Setting.baseDir,[a.Setting.animal '_' a.allAUC.analysisType '.pdf']));

%% AUC R vs NR - plot CDF
layers = {'sup','int','deep'};

for L = 3:length(layers)
close all;

S = a.allAUC.(layers{L}); S.analysisType = a.allAUC.analysisType;
B = a.dBeta.(layers{L}); B.analysisType = a.dBeta.analysisType;

S.layer = a.allAUC.(layers{L}).name;
B.layer = a.dBeta.(layers{L}).name;

[S.h,S.p,S.ks2stat]=kstest2(S.NR,S.R);
[B.h,B.p,B.ks2stat]=kstest2(B.NR,B.R);

%AUC
hAUC = figure(3); hold on;
hAUC_NR_v_R(1) = cdfplot(S.NR);
hAUC_NR_v_R(2) = cdfplot(S.R); 
set(hAUC_NR_v_R,'Linewidth',2)
set(gca,'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02],...
  'xlim',[0 1]);
grid on;
txtSup=['Non-red' newline 'n=' num2str(length(S.NR))]; 
txtInt = ['Red' newline 'n=' num2str(length(S.R))]; 
lAUC = legend(txtSup,...
    txtInt,...
    'Box','off',...
    'Color','none',...
    'Location','Southeast');

strAUC = ['Non-red vs red: p=' num2str(round(S.p,2,'significant'))];
text(0.01,0.9,sprintf(strAUC));

title([S.analysisType ' of NR vs R ' S.layer ' ' a.currentAnimal])
xlabel(S.analysisType);

exportgraphics(gcf,fullfile(a.Setting.baseDir,[a.currentAnimal '_NR_vs_R_'  S.layer '_' S.analysisType '.pdf']));

%Beta
hBeta = figure(4); hold on;
hBeta_NR_v_R(1) = cdfplot(B.NR);
hBeta_NR_v_R(2) = cdfplot(B.R); 
set(hBeta_NR_v_R,'Linewidth',2)
set(gca,'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02]);
grid off;
txtSup=['Non-red' newline 'n=' num2str(length(B.NR))]; 
txtInt = ['Red' newline 'n=' num2str(length(B.R))]; 
lBeta = legend(txtSup,...
    txtInt,...
    'Box','off',...
    'Color','none',...
    'Location','Southeast');

strBeta = ['Non-red vs red: p=' num2str(round(B.p,2,'significant'))];
text(max(xlim),max(ylim)-0.1,sprintf(strBeta),...
    'HorizontalAlignment','Right');

title([B.analysisType ' of NR vs R ' S.layer ' ' a.currentAnimal])
xlabel(B.analysisType);

exportgraphics(gcf,fullfile(a.Setting.baseDir,[a.currentAnimal '_NR_vs_R_' B.layer '_' B.analysisType '.pdf']));

end

%% Beta - plot CDF

[h,a.dBeta.pSI,ks2stat]=kstest2(a.dBeta.sup.NR,a.dBeta.int.NR);
[h,a.dBeta.pSD,ks2stat]=kstest2(a.dBeta.sup.NR,a.dBeta.deep.NR);
[h,a.dBeta.pID,ks2stat]=kstest2(a.dBeta.int.NR,a.dBeta.deep.NR);


hBetaDepth = figure(2); hold on; 
hBeta(1) = cdfplot(a.dBeta.sup.NR);
hBeta(2) = cdfplot(a.dBeta.int.NR); 
hBeta(3) = cdfplot(a.dBeta.deep.NR);
set(hBeta,'Linewidth',2)
set(gca,'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02]);
grid off;
txtSup=['Superficial (< 200 \mum)' newline 'n=' num2str(length(a.dBeta.sup.NR))]; 
txtInt = ['Intermediate (250-300 \mum)' newline 'n=' num2str(length(a.dBeta.int.NR))]; 
txtDeep = ['Deep (>500 \mum)' newline 'n=' num2str(length(a.dBeta.deep.NR))];
lBeta = legend(txtSup,...
    txtInt,...
    txtDeep,...
    'Box','off',...
    'Color','none',...
    'Location','Southeast');

% annotation('textbox',[0.2 0.5 0.3 0.3],'String',['Sup vs Int: ' num2str(pSI) 'Sup vs Deep: ' num2str(pSD) 'Int vs Deep: ' num2str(pID)]);
strBeta = ['Sup vs Int: p=' num2str(round(a.dBeta.pSI,2,'significant')),...
    '\nSup vs Deep: p=' num2str(round(a.dBeta.pSD,2,'significant')),...
    '\nInt vs Deep: p=' num2str(round(a.dBeta.pID,2,'significant'))];
text(max(xlim),max(ylim)-0.1,sprintf(strBeta),...
    'HorizontalAlignment','Right');

title([a.dBeta.analysisType ' of non-' a.Setting.animal ' neurons: comparison across depths'])
xlabel(a.dBeta.analysisType);

exportgraphics(gcf,fullfile(a.Setting.baseDir,[a.Setting.animal '_' a.dBeta.analysisType '.pdf']));

