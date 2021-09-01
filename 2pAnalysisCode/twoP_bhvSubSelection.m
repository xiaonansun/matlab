function bhv = twoP_bhvSubSelection(cBhv)

bhv = cBhv; % Comment out this line when script is converted to function

E = behavior_getStimEvents(bhv);
% [lick, data.lickWinIdx, data.lickWinMs, data.dataLick, data.dataLickTrialNumbers]=twoP_alignToLick(data, SessionData);
[stimTime, stimEndTime, spoutTime, lickR, lickL, levGrabR, levGrabL, water]=behavior_findEventTiming(bhv); 


idxSelfPerformed = ~bhv.SingleSpout & ~bhv.AutoReward;

% Stimulus
bhv.stim.AllRight = idxSelfPerformed & bhv.CorrectSide==2; 
bhv.stim.AllLeft = idxSelfPerformed & bhv.CorrectSide==1;
bhv.stim.EasyRight = idxSelfPerformed & E.ratioEvents' == 1;
bhv.stim.EasyLeft = idxSelfPerformed & E.ratioEvents' == 0;

% Choice
bhv.response.All = idxSelfPerformed & bhv.ResponseSide > 0;
bhv.response.Left = idxSelfPerformed & bhv.ResponseSide == 1;
bhv.response.Right = idxSelfPerformed & bhv.ResponseSide == 2;

% Rewarded trials
bhv.rewarded.All = idxSelfPerformed & bhv.Rewarded; 
bhv.rewarded.Left = idxSelfPerformed & bhv.Rewarded & bhv.ResponseSide == 1;
bhv.rewarded.Right = idxSelfPerformed & bhv.Rewarded & bhv.ResponseSide == 2;

% Error
bhv.error.All = idxSelfPerformed & ~bhv.Rewarded; 
bhv.error.Left = idxSelfPerformed & ~bhv.Rewarded & bhv.ResponseSide == 1;
bhv.error.Right = idxSelfPerformed & ~bhv.Rewarded & bhv.ResponseSide == 2;

bhv.sub.AllIdx = [struct2cell(bhv.stim); struct2cell(bhv.response); struct2cell(bhv.rewarded); struct2cell(bhv.error)];
bhv.sub.AllIdx = vertcat(bhv.sub.AllIdx{:});
bhv.sub.names = [strcat('Stim',fieldnames(bhv.stim));...
    strcat('Response',fieldnames(bhv.response));...
    strcat('Rewarded',fieldnames(bhv.rewarded));...
    strcat('Error',fieldnames(bhv.error))];