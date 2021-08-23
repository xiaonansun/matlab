function rateDisc_testRegs(cPath,Animal,Rec,fileExt,cRegs)
% Code to compute predictive power in different regressor or regressor
% groups.

if ~exist('fileExt','var')
    fileExt = '';
end

if ~exist('cRegs','var') || isempty(cRegs)
    cRegs = []; %test all regressors
end

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end
fPath = [cPath Animal filesep 'SpatialDisc' filesep Rec filesep]; %Widefield data path
load('allenDorsalMapSM.mat', 'dorsalMaps')
allenMask = dorsalMaps.allenMask;
tPath = [fPath 'predVariance' filesep]; %Path to save the results. Make subfolder to get a bit organized.
if ~exist(tPath, 'dir')
    mkdir(tPath);
end

%% load some data and get into right format
betaFile = dir([fPath fileExt 'dimBeta.mat']);
regFile = dir([fPath fileExt 'regData.mat']);
load([fPath betaFile.name],'ridgeVals'); %load model ridge penalty
load([fPath regFile.name],'fullR','regIdx','rejIdx','regLabels') %load design matrix

regIdx = regIdx(~rejIdx);
if sum(fullR(:,1)) == size(fullR,1)
    fullR(:,1) = []; %remove offset if present
    regLabels = regLabels(2:end); %remove offset if present
    regIdx = regIdx(2:end) - 1; %remove offset if present
    rejIdx = rejIdx(2:end); %remove offset if present
end

%get regressor labels and indices for trial segments.
[~, motorLabels, sensorLabels, cogLabels, segIdx] = rateDiscRecordings;

try
    load([fPath 'rsVc.mat'],'U')
catch
    load([fPath 'Vc.mat'],'U')
end
load([fPath 'opts2.mat'],'opts')
opts.frameRate = 15; %make sure this is set to 15Hz
if strcmpi(fileExt, 'org')
    load([fPath 'interpVc.mat'],'Vc','frames') %Vc that was used for the model
else
    load([fPath fileExt 'Vc.mat'],'Vc','frames') %Vc that was used for the model
end

U = alignAllenTransIm(single(U),opts.transParams); %align to allen
U = arrayShrink(U(1:size(allenMask),1:size(allenMask,2),:), allenMask, 'merge');
segFrames = round(segIdx * opts.frameRate);
newVc = rateDisc_getRealignment(fullR, Vc, regIdx, regLabels, segFrames, frames); %re-align Vc to different segments

% assign extra regressor groups that are removed from model together
extraGroups = {'pupils' 'handles' 'licks' 'opMotor' 'Choice' ...
               'handleChoice' 'stimChoice' 'delayChoice' 'respChoice' ...
               'audioStim' 'dlcRegs' 'spontMotor' 'motor' 'sensory' 'cognitive'};
           
extraGroups{2,1} = {'fastPupil' 'slowPupil'};
extraGroups{2,2} = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'};
extraGroups{2,3} = {'lLick' 'rLick'};
extraGroups{2,4} = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'}; %all operant motor regressors
extraGroups{2,5} = {'lhandleChoice' 'rhandleChoice' 'lstimChoice' 'rstimChoice' 'ldelayChoice' 'rdelayChoice' 'lresponseChoice' 'rresponseChoice'}; %all choice regessors

extraGroups{2,6} = {'lhandleChoice' 'rhandleChoice'}; %handle choice regessors
extraGroups{2,7} = {'lstimChoice' 'rstimChoice'}; %stim choice regessors
extraGroups{2,8} = {'ldelayChoice' 'rdelayChoice'}; %delay choice regessors
extraGroups{2,9} = {'lresponseChoice' 'rresponseChoice'}; %response choice regessors

extraGroups{2,10} = {'lfirstAudStim' 'rfirstAudStim' 'lAudStim' 'rAudStim'}; %all auditory stimulus regessors
extraGroups{2,11} = {'cam1_Nosetip' 'cam1_Jaw' 'cam1_Tongue' 'cam2_Nose_tip' 'cam2_Jaw' 'cam2_Tongue' 'cam2_groom'}; %all auditory stimulus regessors
extraGroups{2,12} = motorLabels(~ismember(motorLabels,extraGroups{2,4})); %all spontaneous motor regressors
extraGroups{2,13} = motorLabels; %all motor regressors
extraGroups{2,14} = sensorLabels;
extraGroups{2,15} = cogLabels;

ridgeFolds = 10;    %folds for cross-validation when assessing predicted variance
rng default %reset randum number generator
randIdx = randperm(size(Vc,2)); %generate randum number index if required
% shCurrent: Shuffle current regressor, shOther: Shuffle all other regressors, 
% shOtherMotor: Shuffle all other motor regressors, shOtherSpontMotor: Shuffle all other spontaneous motor regressors, 
% shTaskOtherSpontMotor: Shuffle task and all other spontaneous movements.
% modLabels = {'shCurrent','shOther','shOtherMotor','shOtherSpontMotor','shTaskOtherSpontMotor'}; 
modLabels = {'shCurrent','shOther'}; 
taskLabels = [sensorLabels cogLabels]; %these are all task regressors

%motor regressors, used for shOtherMotor analysis
oMotorLabels = ['pupils' motorLabels(~ismember(motorLabels, extraGroups{2, ismember(extraGroups(1,:),'pupils')}))];
oMotorLabels = ['handles' oMotorLabels(~ismember(oMotorLabels, extraGroups{2, ismember(extraGroups(1,:),'handles')}))];
oMotorLabels = ['licks' oMotorLabels(~ismember(oMotorLabels, extraGroups{2, ismember(extraGroups(1,:),'licks')}))];
oMotorLabels = ['opMotor' oMotorLabels];
oMotorLabels = ['spontMotor' oMotorLabels];

%% test predictive power of individual regressors through cross-validation
Cnt = 0;
if isempty(cRegs)
    testRegs = 0 : length(regLabels) + size(extraGroups,2); %test all regs
else
    testRegs = sort(find(ismember([{'full'} regLabels extraGroups(1,:)],cRegs)))-1; %test selected regs
end

for iRegs = testRegs
    
    Cnt = Cnt+1;
    fprintf('Current regressor is %d of %d\n', Cnt, length(testRegs));

    for modRuns = 1:length(modLabels)
        
        %index for current regressor or group
        if iRegs <= length(regLabels)
            cIdx = regIdx == iRegs; % index for reduced model.
        else
            cIdx = ismember(regIdx, find(ismember(regLabels,extraGroups{2,iRegs - length(regLabels)}))); % index for extra regressor group
        end
        
        %control for shOtherMotor/shOtherSpontMotor condition: Only use for motor regressors.
        checker = true;
        if strcmpi(modLabels{modRuns}, 'shOtherMotor') || strcmpi(modLabels{modRuns}, 'shOtherSpontMotor') || strcmpi(modLabels{modRuns}, 'shTaskOtherSpontMotor')
            if iRegs <= length(regLabels) && iRegs > 0
               checker = any(ismember(regLabels{iRegs},oMotorLabels));
            elseif iRegs > length(regLabels)
               checker = any(ismember(extraGroups{1,iRegs - length(regLabels)},oMotorLabels));
            end
        end
        
        if iRegs == 0
            cLabel = 'full';
        elseif iRegs <= length(regLabels)
            cLabel = regLabels{iRegs};
        else
            cLabel = extraGroups{1,iRegs - length(regLabels)};
        end
        fprintf('Animal: %s, Recording:%s, Regressor:%s\n', Animal, Rec, cLabel)

        if (strcmpi(modLabels{modRuns}, 'shCurrent') || sum(cIdx) > 0) && checker
            
            fakeR = fullR; %copy  design matrix to shuffle up some regressor set
            if strcmpi(modLabels{modRuns}, 'shCurrent') %this is to the shuffle-current regressor
                shIdx = cIdx;
            elseif strcmpi(modLabels{modRuns}, 'shOther') %this is to the shuffle remaining regressors
                shIdx = ~cIdx;
            elseif strcmpi(modLabels{modRuns}, 'shOtherMotor') %this is to shuffle remaining motor regressors
                shIdx = ismember(regIdx, find(ismember(regLabels,taskLabels)));
                shIdx = ~(shIdx | cIdx);
            elseif strcmpi(modLabels{modRuns}, 'shOtherSpontMotor') %this is to the shuffle all spont. motor regressors
                shIdx = ismember(regIdx, find(ismember(regLabels,[taskLabels extraGroups{2,ismember(extraGroups(1,:),{'opMotor'})}]))); %keep task + operant movements
                shIdx = ~(shIdx | cIdx);
            elseif strcmpi(modLabels{modRuns}, 'shTaskOtherSpontMotor') %this is to the shuffle all spont. motor regressors
                shIdx = ismember(regIdx, find(ismember(regLabels,extraGroups{2,ismember(extraGroups(1,:),{'opMotor'})}))); %keep operant movements
                shIdx = ~(shIdx | cIdx);
            end

            %shuffle selected regressors
            for iCol = find(shIdx)
                fakeR(:,iCol) = fullR(randperm(size(fullR,1)),iCol);
            end
            
            %% run cross-validation
            Vm = zeros(size(Vc),'single');
            foldCnt = floor(size(Vc,2) / ridgeFolds);
            
            tic
            for iFolds = 1:ridgeFolds
                tic
                dataIdx = true(1,size(Vc,2));
                
                if ridgeFolds > 1
                    
                    dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
                    
                    [~, betas] = ridgeMML(Vc(:,dataIdx)', fakeR(dataIdx,:), true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                    Vm(:,~dataIdx) = (fakeR(~dataIdx,:) * betas)'; %predict remaining data
                    
                    if rem(iFolds,ridgeFolds/5) == 0
                        fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
                        toc
                    end
                else
                    
                    [~, betas] = ridgeMML(Vc', fakeR, true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                    Vm = (fakeR * betas)'; %predict remaining data
                    disp('Ridgefold is <= 1, fit to complete dataset');
                end
            end
            
            
            %% compute correlation between data and prediction over for all frames
            Vc = reshape(Vc,size(Vc,1),[]);
            Vm = reshape(Vm,size(Vm,1),[]);
            covVc = cov(Vc');  % S x S
            covVm = cov(Vm');  % S x S
            cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
            covP = sum((U * cCovV) .* U, 2)';  % 1 x P
            varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
            varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
            stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
            cMap = (covP ./ stdPxPy)';
            
            %% compute correlation between data and prediction over specific trial segments
            newVm = rateDisc_getRealignment(fullR, Vm, regIdx, regLabels, segFrames, frames); %re-align Vc to different segments
            segMovie = NaN(size(U,1),length(segIdx), 'single');
            for iSegs = 1:length(segIdx)
                try
                    if iSegs == 1
                        cIdx = 1 : sum(segFrames(1:iSegs)); %range of frames for current segment
                    else
                        cIdx = sum(segFrames(1:iSegs-1))+1 : sum(segFrames(1:iSegs)); %range of frames for current segment
                    end
                    
                    % isolate frames for current segment from raw and modeled data
                    cData = newVc(:, cIdx, :);
                    cData = reshape(cData, size(cData,1), []);
                    cData = cData(:, ~isnan(cData(1,:)));
                    
                    mData = newVm(:, cIdx, :);
                    mData = reshape(mData, size(mData,1), []);
                    mData = mData(:, ~isnan(mData(1,:)));
                    
                    covVc = cov(cData');  % S x S
                    covVm = cov(mData');  % S x S
                    cCovV = bsxfun(@minus, mData, mean(mData,2)) * cData' / (size(cData,2) - 1);  % S x S
                    covP = sum((U * cCovV) .* U, 2)';  % 1 x P
                    varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
                    varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
                    stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
                    segMovie(:,iSegs) = (covP ./ stdPxPy)';
                end
            end
            
            %% compute correlation between data and prediction for each frame in all trials
%             cMovie = zeros(size(U,1),frames, 'single');
%             for iFrames = 1:frames
%                 
%                 frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
%                 tic
%                 cData = bsxfun(@minus, Vc(:,frameIdx), mean(Vc(:,frameIdx),2));
%                 cModel = bsxfun(@minus, Vm(:,frameIdx), mean(Vm(:,frameIdx),2));
%                 covVc = cov(cData');  % S x S
%                 covVm = cov(cModel');  % S x S
%                 cCovV = cModel * cData' / (length(frameIdx) - 1);  % S x S
%                 covP = sum((U * cCovV) .* U, 2)';  % 1 x P
%                 varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
%                 varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
%                 stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
%                 cMovie(:,iFrames) = gather(covP ./ stdPxPy)';
%                 clear cData cModel
%                 
%                 if rem(iFrames,round(frames/4)) == 0
%                     fprintf(1, 'Current frame is %d out of %d\n', iFrames,frames);
%                     toc
%                 end
%             end
%             fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));
            
            %% save results           
            if ~exist([tPath modLabels{modRuns}],'dir')
                mkdir([tPath modLabels{modRuns}]);
            end
            save([tPath modLabels{modRuns} filesep fileExt cLabel 'corr.mat'], 'cMap', 'segMovie', 'iRegs', 'regLabels', '-v7.3');
        end
    end
    fprintf('Finished. Current reg: %d\n', Cnt);
end
save([tPath fileExt 'extraGroups.mat'], 'extraGroups'); %save labels for extra groups
save([tPath fileExt 'oMotorLabels.mat'], 'oMotorLabels'); %save other motor labels where some regressors sets are combined