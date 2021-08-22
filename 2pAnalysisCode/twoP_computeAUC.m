function AUC = twoP_computeAUC(AUC,LR)

for i = 1:length(AUC.iSelected)
    
    ii = AUC.iSelected(i);
    iNR = LR(ii).idx_notredcell;
    iR = LR(ii).idx_redcell;
    sAUC = LR(ii).allAUC;
    sMu = LR(ii).shufMu;
    sSigma = LR(ii).shufSigma;
    
    if ~isempty(LR(ii).segFrames)
        segFrames = LR(ii).segFrames;
    else
        segFrames = [30    53    91   106   136];
    end
    
    if ~isempty(LR(ii).cvAcc) 
    AUC.epoch(ii).iPreStim = [segFrames(1) segFrames(2)]; iPS = AUC.epoch(ii).iPreStim;
    AUC.epoch(ii).iEarlyStim = [segFrames(2)+1 floor(mean(segFrames(2)+1:segFrames(3)))]; iES = AUC.epoch(ii).iEarlyStim;
    AUC.epoch(ii).iLateStim = [ceil(mean(segFrames(2)+1:segFrames(3))) segFrames(3)]; iLS = AUC.epoch(ii).iLateStim;
    AUC.epoch(ii).iDelay = [segFrames(3)+1 segFrames(4)]; iD = AUC.epoch(ii).iDelay;
    AUC.epoch(ii).iEarlyRes = [segFrames(1,4)+1 floor(mean(segFrames(4)+1:segFrames(5)))]; iER = AUC.epoch(ii).iEarlyRes;
    AUC.epoch(ii).iLateRes = [ceil(mean(segFrames(4)+1:segFrames(5))) segFrames(5)]; iLR = AUC.epoch(ii).iLateRes;
    AUC.NR.val = [AUC.NR.val; mean(sAUC(iNR,iPS(1):iPS(2)),2,'omitnan') ...
        mean(sAUC(iNR,iES(1):iES(2)),2,'omitnan') ...
        mean(sAUC(iNR,iLS(1):iLS(2)),2,'omitnan') ...
        mean(sAUC(iNR,iD(1):iD(2)),2,'omitnan') ...
        mean(sAUC(iNR,iER(1):iER(2)),2,'omitnan') ...
        mean(sAUC(iNR,iLR(1):iLR(2)),2,'omitnan')];
    AUC.R.val = [AUC.R.val; mean(sAUC(iR,iPS(1):iPS(2)),2,'omitnan') ...
        mean(sAUC(iR,iES(1):iES(2)),2,'omitnan') ...
        mean(sAUC(iR,iLS(1):iLS(2)),2,'omitnan') ...
        mean(sAUC(iR,iD(1):iD(2)),2,'omitnan') ...
        mean(sAUC(iR,iER(1):iER(2)),2,'omitnan') ...
        mean(sAUC(iR,iLR(1):iLR(2)),2,'omitnan')];
    AUC.NR.valAll= [AUC.NR.valAll;sAUC(iNR,:)]; AUC.R.valAll = [AUC.R.valAll;sAUC(iR,:)];
    AUC.NR.mu = [AUC.NR.mu;sMu(iNR,:)]; AUC.R.mu = [AUC.R.mu;sMu(iR,:)];
    AUC.NR.sigma = [AUC.NR.sigma;sSigma(iNR,:)]; AUC.R.sigma = [AUC.R.sigma;sSigma(iR,:)];
    else
        continue
    end
end