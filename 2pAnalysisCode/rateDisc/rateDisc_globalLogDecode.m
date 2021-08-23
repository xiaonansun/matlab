function rateDisc_globalLogDecode(animal, cPath, tPath)
%% check some basic variables
if ~exist('regType','var') || isempty(regType)
    regType = 'lasso';
end

if ~exist('learnType','var') || isempty(learnType)
    learnType = 'leastsquares';
end

if ~exist('stepSize','var') || isempty(stepSize)
    stepSize = 1;
end

if ~exist('useTrials','var') || isempty(useTrials)
    useTrials = 50; %minum nr of trials for decoder
end

cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';

%% 
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
dPath = [cPath animal filesep 'SpatialDisc' filesep]; % path for raw data

%% raw data
load([bPath 'trialInfo.mat'], 'trialCnt', 'recs');
wVfile = matfile([bPath 'wV.mat']); %use this to load temporal dimensions


load([bPath 'wV.mat'],'blockInd'); %use this to rebuild blocks
load([bPath 'wV.mat'],'wU');
load([bPath 'bV.mat'],'bU');
load([bPath 'mask.mat'],'allenMask','xRange','yRange');
allenMask = allenMask(yRange,xRange);

cData = rateDisc_blockRebuild(allenMask,blockInd, bU, wV(1:nDims(iDims),iFrames+1 : iFrames + stepSize), wU(:,1:nDims(iDims)));

%%




    cIdx = sum(trialCnt(1:iRecs-1))+1 : sum(trialCnt(1:iRecs)); %trials for current recording
    wV = wVfile.wV(:,:,cIdx); %load global dimensions
    wV = reshape(wV(:,~isnan(wV(1,:))),size(wV,1),[]); %exclude NaN frames
    wV = bsxfun(@minus, wV, mean(wV,2)); %make sure wV is zero-mean




fPath = [cPath animal filesep 'SpatialDisc' filesep rec filesep]; %Widefield data path
load([fPath 'mask.mat'], 'mask');
load([fPath 'opts.mat'], 'opts');
bhvFile = dir([fPath animal '_SpatialDisc*.mat']);
load([fPath bhvFile(1).name], 'SessionData'); %load behavior data
load([fPath 'Vc.mat'], 'bTrials', 'U'); %spatial dims
load([fPath 'interpVc.mat'], 'Vc', 'frames'); %temporal dims as used in regression model
load([fPath 'interpVmotor.mat'], 'Vmotor'); %reconstruction based on motor regressors
load([fPath 'orgregData.mat'], 'fullR', 'regLabels', 'regIdx', 'rejIdx');

regIdx = regIdx(~rejIdx);
vidR = fullR(:,ismember(regIdx, find(ismember(regLabels,{'Move' 'bhvVideo'}))))'; %video regressors
segFrames = floor(segIdx * opts.frameRate); %max nr of frames per segment
U = arrayShrink(U, mask);

load([fPath 'regData.mat'], 'fullR', 'regLabels', 'regIdx', 'rejIdx', 'spoutR', 'trialIdx');
spoutR(trialIdx,:) = []; %reject bad trials from spoutR

Vc = rateDisc_getRealignment(fullR, spoutR, Vc, regIdx(~rejIdx), regLabels, segFrames, frames); %align to different trial segments
Vmotor = rateDisc_getRealignment(fullR, spoutR, Vmotor, regIdx(~rejIdx), regLabels, segFrames, frames); %align to different trial segments
vidR = rateDisc_getRealignment(fullR, spoutR, vidR, regIdx(~rejIdx), regLabels, segFrames, frames); %align to different trial segments

trialIdx = unique(ceil(find(~trialIdx)/frames)); %change to single trial index (instead of frames)
bhv = selectBehaviorTrials(SessionData,bTrials(trialIdx)); %only use completed trials that are in the Vc dataset

%% run decoder
[cvChoice, cvBeta, trialCnt] = rateDisc_logDecoder(Vc, U, bhv, useTrials, targMod, regType);
[noMotorChoice, noMotorBeta] = rateDisc_logDecoder(Vc - Vmotor, U, bhv, useTrials, targMod, regType);
[motorChoice, motorBeta] = rateDisc_logDecoder(Vmotor, U, bhv, useTrials, targMod, regType);
vidChoice = rateDisc_logDecoder(vidR, [], bhv, useTrials, targMod, regType);

%% save results to local server
tPath = [tPath animal filesep 'SpatialDisc' filesep rec filesep]; %Widefield data path
if ~exist(tPath, 'dir')
    mkdir(tPath);
end

%save results
save([tPath 'logDecode_' regType num2str(targMod) '.mat'], 'trialCnt', 'cvChoice', 'cvBeta', '-v7.3');
save([tPath 'logDecodeMotor_' regType num2str(targMod) '.mat'], 'motorChoice', 'motorBeta', '-v7.3');
save([tPath 'logDecodeNoMotor_' regType num2str(targMod) '.mat'], 'noMotorChoice', 'noMotorBeta', '-v7.3');
save([tPath 'logDecodeVideo_' regType num2str(targMod) '.mat'], 'vidChoice', '-v7.3');

end