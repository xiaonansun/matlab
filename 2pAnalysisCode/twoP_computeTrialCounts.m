function useTrials = twoP_computeTrialCounts(animal,session,useTrials)

%%

% animal = 'CSP30'; session = '200302';

% useTrials = 250;

S = twoP_settings;
imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;

Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mat'),'Vc'); dataIn = Vc.Vc;
cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); bhv = cBhv.cBhv;

perfIdx = ~bhv.DidNotChoose & ~bhv.DidNotLever; %performed trials with target modality

dataIn = dataIn(:, :, perfIdx); %reject non-performed trials

leftIdx = (bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded); %trials were animal went left
pLeftIdx = [false leftIdx(1:end-1)];
pLeftIdx = pLeftIdx(perfIdx); %previous choices
nLeftIdx = [leftIdx(2:end) false];
nLeftIdx = nLeftIdx(perfIdx); %future choices
leftIdx = leftIdx(perfIdx); %current choices

cData = squeeze(dataIn(:, 1, :));

choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, [], useTrials); %equalize L/R choices

selectedTrials = sum(choiceIdx);

while useTrials > selectedTrials
    useTrials = floor(selectedTrials/10)*10;
    choiceIdx = rateDisc_equalizeTrials(~isnan(cData(1,:)), leftIdx, [], useTrials); %equalize L/R choices
    selectedTrials = sum(choiceIdx);
end
    
    
% if selectedTrials < useTrials

disp(['Total number of self-performed trials: ' num2str(sum(perfIdx))...
    '. Total number of selected trials: ' num2str(selectedTrials) '.']);