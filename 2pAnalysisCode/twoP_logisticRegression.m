function lr = twoP_logisticRegression(animal, session, data, SessionData, subSampling, loadExisting)
% animal: name of animal to import
% session: name of session (subfolder) to import
% data: two-photon data struct
% SessionData: behavior data

dbstop if error

baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy'; % directory of all 2p data
lrDir = fullfile(baseDir,animal,'imaging',session,'logisticRegression'); % directory to save the current logistic regression analysis data and figures
if ~exist(lrDir,'dir') % if the current analysis directory doesn't exist, then create one
    mkdir(lrDir);
end
lrDataPath = fullfile(lrDir, [animal '_' session '_lr.mat']); % create a file path for the current session's data (.mat file)

if exist('loadExisting','var') 
    if loadExisting == true
        load(lrDataPath,'lr');
        return
    end
end

% Initialize some parameters for logistic regression
regType = 'lasso'; %lasso or ridge
stepSize = [];
decType = 'allChoice';
segIdx = [1 0.75 1.25 0.5 1];
reps = 10; % number of times to randomly subsample non-red cells
minReps = 10; % number of minimum subsamples for data to be included for analysis
fTrials = 0.5; % fraction of total trials to use for logistic regression
opts.preStim = data.trialStimFrame*data.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
opts.frameRate = 1000/data.msPerFrame; % Frame rate of imaging
sRate = opts.frameRate;
segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment
loadExistingData = false;
cBhv = selectBehaviorTrials(SessionData, data.trialNumbers); %% very important in matching trial indices

lr.cBhv = cBhv;

if exist(lrDataPath,'file') && loadExistingData == true
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

% Aligns imaging data based on the behavior data
Vc = data.neural;
Vc_r = data.neural(data.idx_redcell ,:,:);
Vc_nr = data.neural(data.idx_notredcell,:,:);

% Define the useTrials parameter
useTrials = 150;
% if length(data.trialNumbers) >= 250
%     useTrials = round(size(Vc,3)*fTrials/10)*10;
% end

% if length(data.trialNumbers) >= 500
%     useTrials = 400;
% elseif length(data.trialNumbers) >= 400 && length(data.trialNumbers) < 500
%     useTrials = 300;
% elseif length(data.trialNumbers) >= 300 && length(data.trialNumbers) < 400
%     useTrials = 200;
% elseif length(data.trialNumbers) < 200
%     disp('There are fewer than 200 trials in this session, terminating analysis...');
%     return
% end

Vc = rateDisc_getBhvRealignment(Vc, cBhv, segFrames, opts); % re-aligned imaging data to trial epoches
Vc_r = rateDisc_getBhvRealignment(Vc_r, cBhv, segFrames, opts); 
Vc_nr = rateDisc_getBhvRealignment(Vc_nr, cBhv, segFrames, opts); 


[lr.cvAcc, lr.bMaps, lr.betaNeuron, lr.mdlAll, lr.trialCnt, lr.allAUC, lr.shufAUC] = rateDisc_logDecoder(Vc, [], cBhv, useTrials, 0, regType, stepSize, decType,[],1);
lr.Vc = Vc; lr.Vc_r = Vc_r; lr.Vc_nr = Vc_nr;
while all(isnan(lr.cvAcc))
    fTrials = fTrials-0.1; lr.fTrials = fTrials;
    useTrials = round(size(Vc,3)*fTrials/10)*10;
    lr.useTrials = useTrials; 
    [lr.cvAcc, lr.bMaps, lr.betaNeuron, lr.mdlAll, lr.trialCnt, lr.allAUC] = rateDisc_logDecoder(Vc, [], cBhv, useTrials, 0, regType, stepSize, decType);
end
[lr.cvAcc_r, bMaps_r, betaNeuron_r, mdlAll_r, trialCnt_r] = rateDisc_logDecoder(Vc_r, [], cBhv, useTrials, 0, regType, stepSize, decType); 

%Subsample non-red cells
tic
if subSampling == true
    cvAcc_nr_rep = zeros(reps,size(Vc,2));
    rand_cell_idx = zeros(reps,length(data.idx_redcell));
    for i = 1:reps
        rand_cell_idx(i,:) = randperm(length(data.idx_notredcell),length(data.idx_redcell));
    end
    parfor i = 1:reps
        [cvAcc_nr_rep(i,:), bMaps, betaNeuron, mdlAll, trialCnt_nr_rep(i,:)] = rateDisc_logDecoder(Vc_nr(rand_cell_idx(i,:),:,:), [], cBhv, useTrials, 0, regType, stepSize, decType);
    end
    lr.cvAcc_nr_rep = cvAcc_nr_rep;
    lr.mcvAcc_nr_rep = mean(lr.cvAcc_nr_rep,1);
    
    %Define the error shades
    lr.x = 1:numel(lr.mcvAcc_nr_rep); lr.x2 = [lr.x fliplr(lr.x)];
    lr.uBound = lr.mcvAcc_nr_rep+std(cvAcc_nr_rep); lr.lBound = lr.mcvAcc_nr_rep-std(cvAcc_nr_rep);
    lr.inBetween = [lr.uBound fliplr(lr.lBound)];
end
disp(['Subsampling completed in ' num2str(toc) ' seconds']);

%Define the indices of each epoch
iReal = find(~isnan(lr.cvAcc));
iEpoch(1,:) = [1 iReal(find(diff(iReal)>1)+1)];
iEpoch(2,:) = [iReal(diff(iReal)>1) iReal(end)];
iEpoch(3,:) = [1 find(diff(iReal)>1)+1];
lr.iEpoch = iEpoch;
lr.segFrames = segFrames;

% Remove NaNs for plotting
cvAcc = lr.cvAcc; cvAcc(isnan(cvAcc))=[];
mcvAcc_nr_rep = lr.mcvAcc_nr_rep; mcvAcc_nr_rep(isnan(mcvAcc_nr_rep))=[];
cvAcc_r = lr.cvAcc_r; cvAcc_r(isnan(cvAcc_r))=[];
lr.xm = 1:numel(mcvAcc_nr_rep); xm = lr.xm; x2m = [xm fliplr(xm)]; % defines the x-axis values of the upper and lower bound
uBound = lr.uBound; uBound(isnan(uBound))=[]; 
lBound = lr.lBound; lBound(isnan(lBound))=[];
inBetween = lr.inBetween; inBetween = [uBound fliplr(lBound)];

% Create Plots
figure('Units', 'pixels', ...
    'Position', [100 100 600 250]);
hold on;

hAllNeurons = plot(cvAcc,'k','LineWidth',2); % all neurons
hNonRedNeurons = plot(xm, mcvAcc_nr_rep, 'g', 'LineWidth', 2);
hError = fill(x2m, inBetween, 'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.2);
hRedNeurons = plot(cvAcc_r,'r','LineWidth',2); % only tdT+ neurons

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

title([animal ' ' session ', ' num2str(length(data.trialNumbers)) ' trials']);

% legend([hAllNeurons;hNonRedNeurons;hRedNeurons],['All neurons';'tdT- neurons';'tdT+ neurons'],...
legend({['All neurons, n=' num2str(size(data.neural,1))];['tdT- neurons, n=' num2str(length(data.idx_notredcell))];'';['tdT+ neurons, n=' num2str(length(data.idx_redcell))]},...
    'Box','off',...
    'Color','none',...
    'Location','southeast');

save([lrDir filesep animal '_' session '_lr.mat'],'lr');

exportgraphics(gcf,[lrDir filesep animal '_' session '_lr.pdf']);