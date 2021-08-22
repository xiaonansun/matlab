%% Logistic regression - old code
baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy'; % directory of all 2p data
lrDir = fullfile(baseDir,animal,'imaging',session,'logisticRegression'); % directory to save the current logistic regression analysis data and figures
if ~exist(lrDir,'dir') % if the current analysis directory doesn't exist, then create one
    mkdir(lrDir);
end
lrDataPath = fullfile(lrDir, [animal '_' session '_lr.mat']); % create a file path for the current session's data (.mat file)

% Initialize some parameters for logistic regression
regType = 'lasso'; %lasso or ridge
stepSize = [];
decType = 'allChoice';
segIdx = [1 0.75 1.25 0.5 1];
reps = 40; % number of times to randomly subsample non-red cells
minReps = 10; % number of minimum subsamples for data to be included for analysis
subSampling = true; % Perform subsampling in the current analysis

cBhv = SessionData;

if exist(lrDataPath,'file')
    disp('Existing analyzed data exists, loading...');
    load(lrDataPath)
    if size(lr.cvAcc_nr_rep,1) >= minReps
        subSampling = false;
    else
        disp('There are less than 10 subsamples, will re-run...');
        subSampling = true;
    end
end

if ~exist('subSampling','var') || isempty(decType)
    subSampling = false;
end

opts.preStim = data.trialStimFrame*data.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
opts.frameRate = 1000/data.msPerFrame; % Frame rate of imaging
sRate = opts.frameRate;
segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment

% Aligns imaging data based on the behavior data
Vc = data.neural;
Vc_r = data.neural(data.idx_redcell ,:,:);
Vc_nr = data.neural(data.idx_notredcell,:,:);
% Vc_nr_matched = data.neural(data.idx_notredcell(1:length(data.idx_redcell)),:,:);
% idx_notredcell_shuffled = data.idx_notredcell(randperm(length(data.idx_notredcell)));
% Vc_nr_matched_shuffled = data.neural(idx_notredcell_shuffled(1:length(data.idx_redcell)),:,:);
% Vc_nr_matched_shuffled = data.neural(data.idx_notredcell(randperm(length(data.idx_notredcell),length(data.idx_redcell))),:,:);

useTrials = floor(size(Vc,3)*0.9);
% useTrials = 350;

Vc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts); %aligned to different trial episodes
Vc_r = rateDisc_getBhvRealignment(Vc_r, cBhv, segFrames, opts); %aligned to different trial episodes
Vc_nr = rateDisc_getBhvRealignment(Vc_nr, cBhv, segFrames, opts); %aligned to different trial episodes
% Vc_nr_matched = rateDisc_getBhvRealignment(Vc_nr_matched, cBhv, segFrames, opts); %aligned to different trial episodes
% Vc_nr_matched_shuffled = rateDisc_getBhvRealignment(Vc_nr_matched_shuffled, cBhv, segFrames, opts); %aligned to different trial episodes


% Create plots
% figure('Units', 'pixels', ...
%     'Position', [100 100 500 375]);

% plot(lr.cvAcc,'k','LineWidth',2); % all neurons
% hold on;
% plot(lr.cvAcc_nr,'--g','LineWidth',2.0); 
% plot(lr.cvAcc_r,'r','LineWidth',2); % only tdT+ neurons
% plot(mean(lr.cvAcc_nr_rep,1),'m','LineWidth',2.0)

% The following loop ramdomly selects n = length(data.idx_redcell) neurons
% from non-red cells and repeat i = reps number of times
if subSampling == true
    cvAcc_nr_rep = zeros(reps,size(Vc,2));
    rand_cell_idx = zeros(reps,length(data.idx_redcell));
    for i = 1:reps
        rand_cell_idx(i,:) = randperm(length(data.idx_notredcell),length(data.idx_redcell));
    end
    parfor i = 1:reps
        [cvAcc_nr_rep(i,:), bMaps, trialCnt_nr_rep(i,:)] = rateDisc_logDecoder(Vc_nr(rand_cell_idx(i,:),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType);
    end
    lr.cvAcc_nr_rep = cvAcc_nr_rep;
    lr.mcvAcc_nr_rep = mean(lr.cvAcc_nr_rep);
    lr.x = 1:numel(lr.mcvAcc_nr_rep);
    lr.uBound = lr.mcvAcc_nr_rep+std(cvAcc_nr_rep); lr.lBound = lr.mcvAcc_nr_rep-std(cvAcc_nr_rep);
    lr.x2 = [lr.x fliplr(lr.x)];
    lr.inBetween = [lr.uBound fliplr(lr.lBound)];
%     fill(lr.x2, lr.inBetween, 'g',...
%         'EdgeColor','g');
%     plot(lr.x, lr.mcvAcc_nr_rep, 'g', 'LineWidth', 2);
% else
%     plot(lr.cvAcc_nr_matched_shuffled,'g','LineWidth',2.0); % only tdT- neurons, randomly selected from all tdT- neurons, matched in number to tdT+ neurons
end

[lr.cvAcc, bMaps, trialCnt] = rateDisc_logDecoder(Vc, [], cBhv, useTrials, 0, regType, stepSize, decType); 
[lr.cvAcc_r, bMaps, trialCnt_r] = rateDisc_logDecoder(Vc_r, [], cBhv, useTrials, 0, regType, stepSize, decType); 
[lr.cvAcc_nr, bMaps, trialCnt_nr] = rateDisc_logDecoder(Vc_nr, [], cBhv, useTrials, 0, regType, stepSize, decType); 
% [lr.cvAcc_nr_matched, bMaps, trialCnt_nr_matched] = rateDisc_logDecoder(Vc_nr_matched, [], cBhv, 400, 0, regType, stepSize, decType); 
% [lr.cvAcc_nr_matched_shuffled, bMaps, trialCnt_nr_matched_shuffled] = rateDisc_logDecoder(Vc_nr_matched_shuffled, [], cBhv, 400, 0, regType, stepSize, decType); 

lr.mcvAcc_nr_rep = mean(lr.cvAcc_nr_rep);

iReal = find(~isnan(lr.cvAcc));
iEpoch(1,:) = [1 iReal(find(diff(iReal)>1)+1)];
iEpoch(2,:) = [iReal(diff(iReal)>1) iReal(end)];
iEpoch(3,:) = [1 find(diff(iReal)>1)+1];

% Remove Nan
lr.cvAcc(isnan(lr.cvAcc))=[];
lr.cvAcc_nr(isnan(lr.cvAcc_nr))=[];
% lr.cvAcc_nr_matched(isnan(lr.cvAcc_nr_matched))=[];
% lr.cvAcc_nr_matched_shuffled(isnan(lr.cvAcc_nr_matched_shuffled))=[];
lr.mcvAcc_nr_rep(isnan(lr.mcvAcc_nr_rep))=[];
lr.xm = 1:numel(lr.mcvAcc_nr_rep);
lr.x2m = [lr.xm fliplr(lr.xm)];
lr.cvAcc_r(isnan(lr.cvAcc_r))=[];
lr.uBound(isnan(lr.uBound))=[];
lr.lBound(isnan(lr.lBound))=[];
lr.inBetween = [lr.uBound fliplr(lr.lBound)];

% Create Plots
figure('Units', 'pixels', ...
    'Position', [100 100 500 250]);
hold on;

hAllNeurons = plot(lr.cvAcc,'k','LineWidth',2); % all neurons
hNonRedNeurons = plot(lr.xm, lr.mcvAcc_nr_rep, 'g', 'LineWidth', 2);
hError = fill(lr.x2m, lr.inBetween, 'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.2);
hRedNeurons = plot(lr.cvAcc_r,'r','LineWidth',2); % only tdT+ neurons

labelEpoch = {'pre-stim','stimulus','delay','choice';'top','top','bottom','bottom'};

for i = 1:size(iEpoch,2)
    xline(iEpoch(3,i),'-',labelEpoch{1,i},...
        'LabelVerticalAlignment',labelEpoch{2,i},...
        'FontSize',8);
end

set(gca, ...
  'Box'         , 'off'     , ...
  'Color'       , 'none'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'YTick'       , 0:0.1:1   , ...
  'LineWidth'   , 1);
ylabel('Accuracy'); xlabel('Frame');

legend([hAllNeurons;hNonRedNeurons;hRedNeurons],['All neurons ';'tdT- neurons';'tdT+ neurons'],...
    'Box','off',...
    'Color','none',...
    'Location','east');

save([lrDir filesep animal '_' session '_lr.mat'],'lr');

exportgraphics(gcf,[lrDir filesep animal '_' session '_lr.pdf']);