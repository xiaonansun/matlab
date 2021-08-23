options = struct(); % prepare TUDA options
options.K = 6;
options.DirichletDiag = 100;
options.embeddedlags = -3:3; 
options.initrep = 4; 
options.initcyc = 10; 
options.cyc = 50; 
options.plotAverageGamma = 1;
options.pca = 200;

%%
cd('C:\Users\smusall\Google Drive\mSM63_06-Jul-2018')
load opts2
load Vc;
U = arrayShrink(U,mask,'merge');
load('mSM63_SpatialDisc_Jul06_2018_Session2.mat')
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dat
leftIdx = single((bhv.CorrectSide == 1 & bhv.Rewarded) | (bhv.CorrectSide == 2 & ~bhv.Rewarded)); %trials were animal went left (choice)
leftIdx(leftIdx == 0) = -1;

useIdx = ~isnan(mean(Vc(1,:,:),3));
Vc = Vc(:,useIdx,:);
T = repmat(size(Vc,2),size(Vc,3),1);
% Y = repmat(leftIdx',sum(useIdx),1);
Y = (leftIdx');
X = reshape(Vc,size(Vc,1),[])';

tic; [tuda,Gamma,GammaInit,vpath,stats] = tudatrain(X,Y,T,options); toc
Gamma  = reshape(Gamma, [], length(T), size(Gamma,2));
Betas = squeeze(tudabeta(tuda));

%%
figure; 
subplot(1,2,1); plot(squeeze(nanmean(Gamma,2)),'linewidth',2); ylim([0 1]); xlim([4 size(Gamma,1)-3]); axis square
subplot(1,2,2); plot(squeeze((stats.R2_states)),'linewidth',2); ylim([-0.2 1]); xlim([1 size(Gamma,1)]); axis square

%% sparsify weights
Gamma1  = reshape(Gamma, [], size(Gamma,3));
[tuda,Gamma1,encmodel,decmodel] = tudasparsify(X,Y,T,tuda,Gamma1,0.2,0.5);

%% get events for each trial
taskEvents = false(size(Vc,2),size(Vc,3),4);
taskEventType = repmat(3,1,4);

for iTrials = 1 : length(trials)
    stimOn = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    handleOn = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    clear cIdx
    cIdx(1) = handleOn(find(handleOn == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    cIdx(2) = stimOn;
    cIdx(3) = max(cat(2,bhv.stimEvents{iTrials}{:})) + stimOn; %time of last stimulus event
    cIdx(4) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1);
    cIdx = floor((cIdx - stimOn + opts.preStim) * opts.frameRate); %convert to frames. This is the last frame of each segment.
    for x = 1 : length(cIdx)
        taskEvents(cIdx(x),iTrials,x) = true;
    end
end

%% asign some basic options for the model
opts.mPreTime = ceil(2 * opts.frameRate);  % precede motor events to capture preparatory activity in frames
opts.mPostTime = ceil(2 * opts.frameRate);   % follow motor events for mPostStim in frames
opts.framesPerTrial = size(Vc,2); % nr. of frames per trial

Gamma1 = padarray(Gamma,[3 0 0]);
Gamma1  = reshape(Gamma1, [], size(Gamma,3));

[taskR, taskIdx] = makeDesignMatrix(taskEvents, taskEventType, opts); %make design matrix for task variables
[ridgeVals, dimBeta] = ridgeMML(Gamma1, taskR, true); %get ridge penalties and beta weights.
% linBeta = dimBeta(2:end,:); %don't use intercept
% linBeta = dimBeta + mean(reshape(Gamma,[],size(Gamma,3))); %don't use intercept
linBeta = dimBeta; %don't use intercept

%% show response kernels
figure;
for y = 1 : length(unique(taskIdx))
    subplot(2,ceil(length(unique(taskIdx))/2),y);
    hold on
    for x = 1 : size(Gamma,3)
        plot(linBeta(taskIdx == y,x), 'linewidth',2); axis square
    end
end

%%
    
% Vc = rateDisc_getBhvRealignment(Vc, bhv, segFrames, opts); %aligned to different trial episodes

uB = U * squeeze(Betas);
compareMovie(arrayShrink(uB,mask,'split'));

%% figure
% uB = U * squeeze(Betas);
% for x = 1 : size(Betas,2)
%     subplot(3, ceil(size(Betas,2) / 3), x);
%     imagesc(arrayShrink(uB(:,x),mask,'split')); axis image;
% end
%


%%
tic; [acc,acc_star,Ypred,Ypred_star] = tudacv(X,Y,T,options); toc

