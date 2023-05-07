%% Loads the outputs of the encoding model for further analysis/visualization
clear all;
animal = 'CSP29'; session = '20200305a';
% animal = 'Plex50'; session = '200401a';
lFile={'interpVc'; 'orgdimBeta';'orgfullcorr';'orgregData'};

fileExt='.mat';
S = twoP_settings;
dataDir = fullfile(S.dir.imagingRootDir,animal,'imaging',session,'encodingModel');

if ~exist(dataDir,'dir')
dataRootDir='H:\twoP';
dataDir = fullfile(dataRootDir,animal,'SpatialDisc',session);
end

for i = 1:length(lFile)
    %% i=1
    load(fullfile(dataDir,[lFile{i} fileExt]));
end


%%
rho=1;
filterFrames=5;
numOfTaskRegs = 26;
% regNameIdx=25;
betaKernel=struct;
taskIdx = [ones(1,numOfTaskRegs) zeros(1,length(regLabels)-numOfTaskRegs)];
taskBeta = dimBeta(ismember(regIdx,find(taskIdx)),:);
for regNameIdx = 1:length(regLabels)
    display(['Regressor: ' regLabels{regNameIdx}])
    try
        betaKernel(regNameIdx).value=dimBeta(regIdx(~rejIdx) == regNameIdx,:)';
        betaKernel(regNameIdx).regressorName=regLabels{regNameIdx};
    catch ME
        disp(ME.message)
    end
end
betaIdx = find(arrayfun(@(betaKernel) ~isempty(betaKernel.value),betaKernel));

% create regressor summary table
regSummaryTable = table({betaKernel.regressorName}',cellfun(@(x) size(x,2),{betaKernel.value})');
regSummaryTable.Properties.VariableNames = {'RegressorName';'NumOfFrames'};
%% Unsorted beta kernels

figure(1);
hTile = tiledlayout(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))));
hTile.TileSpacing = 'tight';
colors = cbrewer('seq','Greys',1024);

% subtitle(animal)
for i = 1:length(betaIdx)
    nexttile
    regBeta=betaKernel(betaIdx(i)).value;
    imagesc(regBeta,[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    if sum(i == 1:hTile.GridSize(2):(hTile.GridSize(1)*hTile.GridSize(2))) == 0; set(gca,'YTickLabel',[]); end % remove y tick labels on all except the first column
    colormap(colors);
    title(regLabels{betaIdx(i)});
end

%% Sorted beta kernels
figure(2);
hTile = tiledlayout(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))));
hTile.TileSpacing = 'tight';
colors = cbrewer('seq','Greys',1024);
sPerFrame = S.msPerFrame/1000;

for i = 1:length(betaIdx)
    nexttile
    regBeta=betaKernel(betaIdx(i)).value;
%     time_vec = 0:sPerFrame:sPerFrame*size(regBeta,2); % time vector
    [~,iM]=max(smoothdata(regBeta,2,'sgolay',filterFrames),[],2);
    iM=sortrows([iM (1:1:length(iM))'],1);
    imagesc(regBeta(iM(:,2),:),[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    if sum(i == 1:hTile.GridSize(2):(hTile.GridSize(1)*hTile.GridSize(2))) == 0; set(gca,'YTickLabel',[]); end % remove y tick labels on all except the first column
    colormap(colors);
    title(regLabels{betaIdx(i)});
    fig_configAxis(gca);
end

%% Single neuron beta kernels 
id_neuron = 63;
s_per_frame = S.msPerFrame/1000;
figure(3);
subtitle(animal)
y_limits=[min(taskBeta(:,id_neuron)) max(taskBeta(:,id_neuron))];
tileFigure = tiledlayout(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))));
for i = 1:length(betaIdx)
    nexttile
    regBeta=betaKernel(betaIdx(i)).value(id_neuron,:);
    time_vec = 0:s_per_frame:s_per_frame*(length(regBeta)-1);
    line(time_vec,regBeta);
%     plot(regBeta);
    ylim([y_limits(1) y_limits(2)]);
    title(regLabels{betaIdx(i)});
    fig_configAxis(gca);
    offsetAxes(gca);
end
title(tileFigure, ['Neuron #' num2str(id_neuron)])


%% Compute cvR^2 of individual neurons
frames_per_trial = 155;
num_of_trials = size(Vc,2)/frames_per_trial;
neuRsq = single(zeros(size(Vc,1),1));
neu_tr_Rsq = single(zeros(size(Vc,1),num_of_trials));

for i = 1:size(Vc,1)
    cc = corrcoef(Vc(i,:),Vm(i,:)); % compute the correlation coefficient between the model versus data
    neuRsq(i) = cc(1,2)^2; % compute the R^2 for each neuron across all sessions
    for j = 1:num_of_trials
        trial_frames = frames_per_trial*(j-1)+1:frames_per_trial*j;
        neu_tr_Vc = Vc(i,trial_frames);
        neu_tr_Vm = Vm(i,trial_frames);
        cc = corrcoef(neu_tr_Vc,neu_tr_Vm); % compute the single-neuron, single-trial correlation coefficient between the model versus data
        neu_tr_Rsq(i,j) = cc(1,2)^2;
    end
end
neuron_rank = repmat((1:size(neu_tr_Rsq,1))',size(neu_tr_Rsq,2),1);
trial_rank = repmat(1:size(neu_tr_Rsq,2),size(neu_tr_Rsq,1),1); trial_rank = trial_rank(:);
[sorted_neu_tr_Rsq,idx_sorted_neu_tr_Rsq]=sort(neu_tr_Rsq(:),'descend');
sorted_neuron_rank = neuron_rank(idx_sorted_neu_tr_Rsq);
sorted_trial_rank = trial_rank(idx_sorted_neu_tr_Rsq);
sorted_neu_tr_id = [sorted_neuron_rank sorted_trial_rank];

[rankRsq,idxRsq] = sort(neuRsq,'descend');
plot(rankRsq);
fig_configAxis(gca);
offsetAxes(gca);
xlabel('Neurons'); ylabel('cvR^2'); title('Full Model');

figure
imagesc(neu_tr_Rsq);
colors = cbrewer('seq','Greys',1024);
colormap(colors);
xlabel('Trials'); ylabel('Neurons'); title('Single-neuron, single-trial cvR^2')

%% Plot single neuron responses
ranked_id = 501;
id_neuron = sorted_neu_tr_id(ranked_id,1);
id_trial = sorted_neu_tr_id(ranked_id,2);
s_per_frame = S.msPerFrame/1000;
frames_per_trial = 155;
trial_frames = frames_per_trial*(id_trial-1)+1:frames_per_trial*id_trial;
time_vec = 0:s_per_frame:s_per_frame*(length(trial_frames)-1);
neu_tr_Vc = Vc(id_neuron,trial_frames);
neu_tr_Vm = Vm(id_neuron,trial_frames);
figure; hold on;
line(time_vec,neu_tr_Vc);
line(time_vec,neu_tr_Vm);
title(['cvR^2 = ' num2str(neu_tr_Rsq(id_neuron,id_trial)) ...
    ' of neuron #' num2str(id_neuron) ...
    ', trial #' num2str(id_trial)])

%% View the single-trial design matrix of one regressor
% to view regressor names, view regSummaryTable
figure
trial_index = 130;
reg_name = 'rAudStim';
reg_index = find(contains(regLabels,reg_name));
design_matrix_cols = sum(regSummaryTable.NumOfFrames(1:reg_index-1))+1:sum(regSummaryTable.NumOfFrames(1:reg_index));
design_matrix_rows = frames*(trial_index-1)+1:frames*trial_index;
imagesc(fullR(design_matrix_rows,design_matrix_cols));
xlabel('Regressor frames'); ylabel('Trial frames');