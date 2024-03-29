function AUC = twoP_computeAUCFraction(LR,AUC)

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
AUC.NR.pos = AUC.NR.mu + 2*AUC.NR.sigma; AUC.R.pos = AUC.R.mu + 2*AUC.R.sigma;
AUC.NR.neg = AUC.NR.mu - 2*AUC.NR.sigma; AUC.R.neg = AUC.R.mu - 2*AUC.R.sigma; 

% Compute the fraction of choice-selective neurons 
AUC.NR.posSel = AUC.NR.valAll > AUC.NR.pos; AUC.R.posSel = AUC.R.valAll > AUC.R.pos;
AUC.NR.negSel = AUC.NR.valAll < AUC.NR.neg; AUC.R.negSel = AUC.R.valAll < AUC.R.neg;
AUC.R.frPos = sum(AUC.R.posSel,1)/size(AUC.R.posSel,1);
AUC.NR.frPos = sum(AUC.NR.posSel,1)/size(AUC.NR.posSel,1);
AUC.R.frNeg = sum(AUC.R.negSel,1)/size(AUC.R.negSel,1);
AUC.NR.frNeg = sum(AUC.NR.negSel,1)/size(AUC.NR.negSel,1);

AUC.NR.posVal = AUC.NR.valAll.*(AUC.NR.valAll > AUC.NR.pos); AUC.R.posVal = AUC.R.valAll.*(AUC.R.valAll > AUC.R.pos);
AUC.NR.negVal = AUC.NR.valAll.*(AUC.NR.valAll < AUC.NR.neg); AUC.R.negVal = AUC.R.valAll.*(AUC.R.valAll < AUC.R.neg);
AUC.NR.meanPos = sum(AUC.NR.posVal,1,'omitnan')./sum(AUC.NR.posVal~=0,1,'omitnan'); AUC.R.meanPos = sum(AUC.R.posVal,1,'omitnan')./sum(AUC.R.posVal~=0,1,'omitnan');
AUC.NR.meanNeg = sum(AUC.NR.negVal,1,'omitnan')./sum(AUC.NR.negVal~=0,1,'omitnan'); AUC.R.meanNeg = sum(AUC.R.negVal,1,'omitnan')./sum(AUC.R.negVal~=0,1,'omitnan');