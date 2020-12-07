thresholdFound=0;
PM=PAL_AMPM_setupPM('priorAlphaRange',[0:0.01:1],'priorBetaRange',betas,'priorGammaRange',0.5,'priorLambdaRange',0.01,'stimRange',[0:0.05:1],'marginalize','lapse','marginalize','guess','marginalize','slope','numTrials',10000);
mu=0.64;

while(~thresholdFound)
PM=PAL_AMPM_updatePM(PM,detect(mu,PM.xCurrent));
if length(PM.threshold)>10
if abs(PM.threshold(end)-mu)<0.01
thresholdFound=1;
PM.threshold(end)
length(PM.threshold)
end
end
end