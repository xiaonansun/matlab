function sBhv = twoP_bhvSelectTrials(bhv)
%%
idx = ~bhv.SingleSpout & ~bhv.AutoReward & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.StimType == 2; %only use active audio trials
idxR = (bhv.CorrectSide == 1 & bhv.Punished) | (bhv.CorrectSide == 2 & bhv.Rewarded); %right-choice trials
rIdx = rIdx(idx);
