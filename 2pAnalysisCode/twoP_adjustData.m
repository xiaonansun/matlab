function data = twoP_adjustData(data,bhv)

% Sometimes the stimulus-aligned 2P data array contains more trials (3rd
% dimension) than the number of trial codes received. This code discards those extra trials, 
% which is usually one extra trial at the end of the session.
% THIS CODE ONLY NEEDS THE PROCESSED 2P IMAGING STRUCTURE (data)
if length(data.trialNumbers) < size(data.neural,3) 
    disp(['Note: there are ' num2str(size(data.neural,3) - length(data.trialNumbers)) ' fewer data.trialNumbers than size(neural.data,3).'])
    delIdx = find(~ismember(1:1:size(data.neural,3), 1:1:length(data.trialNumbers)));
    data.neural(:,:,delIdx)=[];
    data.stimFramesOrig(delIdx)=[];
    data.DS(:,:,delIdx)=[];
    data.stimSamplesOrig(delIdx)=[];
    data.analog(:,:,delIdx)=[];
    disp(['There were ' num2str(length(delIdx)) ' additional trials in data.neural than in data.trialNumbers']);
else
    disp('The number of trials in data.neural is equal to the number of trials in data.trialNumbers.');
end

% Compares the number of imaging trials to the the number of Bpod trials and discard extra imaging trials
% The number of trials logged by MScan (through the analog trial codes) may
% exceed the number of rewarded trials recorded by Bpod
if max(data.trialNumbers) > length(bhv.Rewarded) 
    disp(['Note: there are ' num2str(max(data.trialNumbers) - length(bhv.Rewarded)) ' fewer rewarded trial(s) than max(data.trialNumbers).']);
    data.trialNumbersOrig = data.trialNumbers;
    delIdx = find(~ismember(data.trialNumbers, 1:1:length(bhv.Rewarded)));
    data.trialNumbers(delIdx)=[];
    data.neural(:,:,delIdx)=[];
    data.stimFramesOrig(delIdx)=[];
    data.DS(:,:,delIdx)=[];
    data.stimSamplesOrig(delIdx)=[];
    data.analog(:,:,delIdx)=[];
else
    disp('The total number of trials captured by MScan equals those recorded by Bpod');
end