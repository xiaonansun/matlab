function rateDisc_logRegress(cPath, tPath, animal, rec, regType, targMod, useTrials)

%% check basic variables
if ~exist('regType','var') || isempty(regType)
    regType = 'lasso';
end

if ~exist('targMod','var') || isempty(targMod)
    targMod = 0;
end

if ~exist('useTrials','var') || isempty(useTrials)
    useTrials = 250; %minum nr of trials for decoder
end
[~, ~, ~, ~, segIdx] = rateDiscRecordings;
stepSize= 3;

%% load some data
fPath = [cPath animal filesep 'SpatialDisc' filesep rec filesep]; %Widefield data path
load([fPath 'opts2.mat'], 'opts');
bhvFile = dir([fPath animal '_SpatialDisc*.mat']);
load([fPath bhvFile(1).name], 'SessionData'); %load behavior data
try
    load([fPath 'rsVc.mat'], 'U'); %spatial dims
catch
    load([fPath 'Vc.mat'], 'U'); %spatial dims
end
mask = isnan(U(:,:,1));
U = arrayShrink(U, mask, 'merge');
opts.frameRate = 15;
load([fPath 'interpVc.mat'], 'bTrials', 'Vc', 'frames'); %get weights of motor model
load([fPath 'vidRegData.mat'], 'vidR'); % get vide regressors
load([fPath 'spontMotorBeta.mat'], 'spontMotorBeta'); %get weights of motor model
load([fPath 'spontMotorregData.mat'], 'spontMotorR'); % get motor design matrix
load([fPath 'motorBeta.mat'], 'motorBeta'); %get weights of motor model
load([fPath 'motorregData.mat'], 'motorR'); % get motor design matrix
spontMotorBeta = nanmean(cat(3,spontMotorBeta{:}),3);
VspontMotor = (spontMotorR * spontMotorBeta)';

motorBeta = nanmean(cat(3,motorBeta{:}),3);
Vmotor = (motorR * motorBeta)';

load([fPath 'orgregData.mat'], 'fullR', 'regLabels', 'regIdx', 'rejIdx','trialIdx');
regIdx = regIdx(~rejIdx);
segFrames = floor(segIdx * opts.frameRate); %max nr of frames per segment

% re-align data
Vc = rateDisc_getRealignment(fullR, Vc, regIdx, regLabels, segFrames, frames); %align to different trial segments
Vmotor = rateDisc_getRealignment(fullR, Vmotor, regIdx, regLabels, segFrames, frames); %align to different trial segments
VspontMotor = rateDisc_getRealignment(fullR, VspontMotor, regIdx, regLabels, segFrames, frames); %align to different trial segments
vidR = rateDisc_getRealignment(fullR, vidR', regIdx, regLabels, segFrames, frames); %align to different trial segments

trialIdx = unique(ceil(find(~trialIdx)/frames)); %change to single trial index (instead of frames)
bhv = selectBehaviorTrials(SessionData,bTrials(trialIdx)); %only use completed trials that are in the Vc dataset

%% run decoder
tic
[cvACC(:,1), cvBeta(:,:,1), trialCnt(:,1)] = rateDisc_logDecoder(Vc, U, bhv, useTrials, targMod, regType, stepSize, 'allChoice'); %all data
[nmACC(:,1), nmBeta(:,:,1)] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, useTrials, targMod, regType, stepSize, 'allChoice'); %no motor
[nsmACC(:,1), nsmBeta(:,:,1)] = rateDisc_logDecoder(Vc - VspontMotor, U, bhv, useTrials, targMod, regType, stepSize, 'allChoice'); %no spont motor

[cvACC(:,2), cvBeta(:,:,2)] = rateDisc_logDecoder(Vc, U, bhv, useTrials, targMod, regType, stepSize, 'allStim'); %all data
[nmACC(:,2), nmBeta(:,:,2)] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, useTrials, targMod, regType, stepSize, 'allStim'); %no motor
[nsmACC(:,2), nsmBeta(:,:,2)] = rateDisc_logDecoder(Vc - VspontMotor, U, bhv, useTrials, targMod, regType, stepSize, 'allStim'); %no spont motor

[cvACC(:,3), cvBeta(:,:,3), trialCnt(:,2)] = rateDisc_logDecoder(Vc, U, bhv, 150, targMod, regType, stepSize, 'Choice'); %all data
[nmACC(:,3), nmBeta(:,:,3)] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, 150, targMod, regType, stepSize, 'Choice'); %no motor
[nsmACC(:,3), nsmBeta(:,:,3)] = rateDisc_logDecoder(Vc - VspontMotor, U, bhv, 150, targMod, regType, stepSize, 'Choice'); %no spont motor

[cvACC(:,4), cvBeta(:,:,4)] = rateDisc_logDecoder(Vc, U, bhv, 150, targMod, regType, stepSize, 'stim'); %all data
[nmACC(:,4), nmBeta(:,:,4)] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, 150, targMod, regType, stepSize, 'stim'); %no motor
[nsmACC(:,4), nsmBeta(:,:,4)] = rateDisc_logDecoder(Vc - VspontMotor, U, bhv, 150, targMod, regType, stepSize, 'stim'); %no spont motor

[cvACC(:,5), cvBeta(:,:,5)] = rateDisc_logDecoder(Vc, U, bhv, 150, targMod, regType, stepSize, 'preChoice'); %all data
[nmACC(:,5), nmBeta(:,:,5)] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, 150, targMod, regType, stepSize, 'preChoice'); %no motor
[nsmACC(:,5), nsmBeta(:,:,5)] = rateDisc_logDecoder(Vc - VspontMotor, U, bhv, 150, targMod, regType, stepSize, 'preChoice'); %no spont motor

[cvACC(:,6), cvBeta(:,:,6)] = rateDisc_logDecoder(Vc, U, bhv, 150, targMod, regType, stepSize, 'nextChoice'); %all data
[nmACC(:,6), nmBeta(:,:,6)] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, 150, targMod, regType, stepSize, 'nextChoice'); %no motor
[nsmACC(:,6), nsmBeta(:,:,6)] = rateDisc_logDecoder(Vc - VspontMotor, U, bhv, 150, targMod, regType, stepSize, 'nextChoice'); %no spont motor

vidChoice(:,1) = rateDisc_logDecoder(vidR, [], bhv, useTrials, targMod, regType, stepSize, 'allChoice');
vidChoice(:,2) = rateDisc_logDecoder(vidR, [], bhv, useTrials, targMod, regType, stepSize, 'stim');
vidChoice(:,3) = rateDisc_logDecoder(vidR, [], bhv, useTrials, targMod, regType, stepSize, 'preChoice');
vidChoice(:,4) = rateDisc_logDecoder(vidR, [], bhv, useTrials, targMod, regType, stepSize, 'nextChoice');

%% save results to local server
tfPath = [tPath animal filesep 'SpatialDisc' filesep rec filesep]; %Widefield data path

%save results
save([tfPath 'logDecode_' regType '_' num2str(stepSize) '_' num2str(targMod) '.mat'], 'trialCnt', 'cvACC', 'cvBeta', '-v7.3');
save([tfPath 'logNoMotorDecode_' regType '_' num2str(stepSize) '_' num2str(targMod) '.mat'], 'trialCnt', 'nmACC', 'nmBeta', '-v7.3');
save([tfPath 'logNoSpontMotorDecode_' regType '_' num2str(stepSize) '_' num2str(targMod) '.mat'], 'trialCnt', 'nsmACC', 'nsmBeta', '-v7.3');
save([tfPath 'logVideoDecode_' regType '_' num2str(stepSize) '_' num2str(targMod) '.mat'], 'vidChoice', '-v7.3');

end