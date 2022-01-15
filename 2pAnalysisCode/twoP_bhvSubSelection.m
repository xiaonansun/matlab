function bhv = twoP_bhvSubSelection(cBhv)
%%
bhv = cBhv; % Comment out this line when script is converted to function

E = behavior_getStimEvents(bhv);
% [lick, data.lickWinIdx, data.lickWinMs, data.dataLick, data.dataLickTrialNumbers]=twoP_alignToLick(data, SessionData);
[stimTime, stimEndTime, spoutTime, lickR, lickL, levGrabR, levGrabL, water]=behavior_findEventTiming(bhv); 


% idxSelfPerformed = ~bhv.SingleSpout & ~bhv.AutoReward;
idxSelfPerformed = ~bhv.SingleSpout & ~bhv.AutoReward & ~bhv.DidNotChoose & ~bhv.DidNotLever & logical(bhv.Assisted) & bhv.StimType == 2; %only use active audio trials

% Stimulus
bhv.stim.AllLeft = idxSelfPerformed & bhv.CorrectSide==1;
bhv.stim.AllRight = idxSelfPerformed & bhv.CorrectSide==2; 
bhv.stim.EasyLeft = idxSelfPerformed & E.ratioEvents' == 0;
bhv.stim.EasyRight = idxSelfPerformed & E.ratioEvents' == 1;

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

bhv.sub.names.stim = strcat('Stim',fieldnames(bhv.stim));
bhv.sub.names.response = strcat('Response',fieldnames(bhv.response));
bhv.sub.names.rewarded = strcat('Rewarded',fieldnames(bhv.rewarded));
bhv.sub.names.error = strcat('Error',fieldnames(bhv.error));
fn = fieldnames(bhv.sub.names);

bhv.sub.AllNames = [];
for i = 1:length(fn)
bhv.sub.AllNames = [bhv.sub.AllNames;bhv.sub.names.(fn{i})];
end
