function [lick,lickWinIdx,lickWinMs,dataLick, dataLickTrialNumbers]=twoP_alignToLick(data, bhv)
% This code uses the stimulus aligned two-photon data matrix (data.neural) and performs
% realignment to the onset of the first lick. 

% Sets some constants
fps=1000/data.msPerFrame; % frames per second
stdThresh=2.5; % Reject delayed lick responses above a threshold, expressed in units of standard deviations from the mean
minLicks=2; % minimum number of licks required to qualify as a response
stimOnsetIdx=find(data.neuralTimes==0); % frame index on which stimulus initiates

%Use the findEventTiming function to extract the timing of stimulus onset, spout movement, and licks.
[stimTime, stimEndTime, spoutTime, lickR, lickL, levGrabR, levGrabL, water]=behavior_findEventTiming(bhv); 
stimTime=stimTime(data.trialNumbers); % Timing of stimulus onset relative to handle grab

% if bhv.nTrials > length(water); water=[water zeros(1,bhv.nTrials-length(water))]; end
% water=water(data.trialNumbers);

if length(bhv.Rewarded) > length(lickR); lickR = [lickR num2cell(zeros(1,bhv.nTrials-length(lickR)))]; end % repairs trial count mismatch
if length(bhv.Rewarded) > length(lickL); lickL = [lickL num2cell(zeros(1,bhv.nTrials-length(lickL)))]; end % repairs trial count mismatch
lickR=lickR(data.trialNumbers); lickL=lickL(data.trialNumbers); % only include trials with imaging data that can be analyzed

lick= zeros(3,length(data.trialNumbers));

% Generate a row of response side: 1 is left and 2 is right
lick(1,:) = bhv.ResponseSide(data.trialNumbers);

% Generate a row of the timing (in seconds) of the first lick event
% relative to trial onset
leftLickIdx=find(cellfun(@numel,lickL)>=minLicks); % Finds the trial 
rightLickIdx=find(cellfun(@numel,lickR)>=minLicks);
[val, ir, il]=intersect(rightLickIdx,leftLickIdx);
rightLickIdx(ir)=[]; leftLickIdx(il)=[]; % remove trials where the animal licked both spouts more than once
leftLickTiming=cellfun(@(x)x(1),lickL((leftLickIdx)));
rightLickTiming=cellfun(@(x)x(1),lickR(rightLickIdx)); %compute the timing of left or right licks
lick(2,leftLickIdx)=leftLickTiming; lick(2,rightLickIdx)=rightLickTiming; % compute the timing of the first lick of a true lick response
lick(2,lick(2,:)==0)=NaN;
lick(2,lick(2,:)>nanmean(lick(2,:))+nanstd(lick(2,:))*stdThresh)=NaN;

% compute the timing of the first lick of a true lick response relative to
% the onset of the sensory stimulus
lick(3,:)=lick(2,:)-stimTime;

% Compute the time index of the first lick
lick(4,:)= round(lick(3,:)*fps)+stimOnsetIdx;
maxPostLickIdx=max(length(data.neuralTimes)-lick(4,:));

lick(5,:)= data.trialNumbers; % Trial IDs; 
lick(6,:)=bhv.Rewarded(lick(5,:)); % Rewarded sessions
lick(7,:)=~bhv.Rewarded(lick(5,:)); % Non-rewarded sessions
lick(8,:)=bhv.AutoReward(lick(5,:)); % Auto-reward sessions

rewardedTrials = lick(5,logical(lick(6,:))); % Index of rewarded trials
lick(:,isnan(lick(1,:)))=[]; % Drops trials  where a NaN was reigstered in the ResponseSide field
[row,col,v] = find(isnan(lick)); lick(:,col)=[]; % Similar to the previous line, this removes trials that contains NaN

preLickFrames = round(fps); postLickFrames = size(data.neural,2)-max(lick(4,:)); % 1 second windows before stimulus
lickWinIdx = -preLickFrames:1:postLickFrames; % IMPORTANT: frame index relative to animal lick
lickWinMs = lickWinIdx*data.msPerFrame;  % IMPORTANT: time (in milliseconds) relative to animal lick
startIdx = lick(4,:)'- preLickFrames -1 ;  endIdx = lick(4,:)' + postLickFrames;
dataLick = zeros(size(data.neural,1),endIdx(1)-startIdx(1),length(lick));
[dataLickTrialNumbers, iL, iD]=intersect(lick(5,:),data.trialNumbers); % <---- need to include iD as an output!!!
for i = 1:length(lick); dataLick(:,:,i) = data.neural(:,startIdx(i):endIdx(i)-1,iD(i)); end
