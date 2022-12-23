%% Loads the outputs of the encoding model for further analysis/visualization

clear all;

animal = 'CSP27'; session = '20200307a';
% session='28-Mar-2020_1';

S = twoP_settings;
dataDir = fullfile(S.dir.imagingRootDir,animal,'imaging',session,'encodingModel');

if ~exist(dataDir,'dir')
dataRootDir='H:\twoP';
dataDir = fullfile(dataRootDir,animal,'SpatialDisc',session);
end

lFile={'interpVc'; 'orgdimBeta';'orgfullcorr';'orgregData'};
fileExt='.mat';

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
% figure(1);
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
subtitle(animal)
for i = 1:length(betaIdx)
    regBeta=betaKernel(betaIdx(i)).value;
    [~,iM]=max(smoothdata(regBeta,2,'sgolay',filterFrames),[],2);
    % [M iM]=max(betaKernel,[],2);
    iM=sortrows([iM (1:1:length(iM))'],1);
    subplot(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))),i)
    imagesc(regBeta,[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    title(regLabels{betaIdx(i)});
end

%% Sorted beta kernels
figure(2);
% suptitle(animal)tiledlayout(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))));
tiledlayout(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))));

colors = cbrewer('seq','Greys',1024);

for i = 1:length(betaIdx)
    nexttile
    regBeta=betaKernel(betaIdx(i)).value;
    [~,iM]=max(smoothdata(regBeta,2,'sgolay',filterFrames),[],2);
    % [M iM]=max(betaKernel,[],2);
    iM=sortrows([iM (1:1:length(iM))'],1);
    imagesc(regBeta(iM(:,2),:),[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    colormap(colors);
    title(regLabels{betaIdx(i)});
end

%% Sorted beta kernels
id_neuron = 138;

figure(3);
subtitle(animal)
y_limits=[min(taskBeta(:,id_neuron)) max(taskBeta(:,id_neuron))];
tiledlayout(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))));
for i = 1:length(betaIdx)
    nexttile
    regBeta=betaKernel(betaIdx(i)).value(id_neuron,:);
%     [~,iM]=max(smoothdata(regBeta,2,'sgolay',filterFrames),[],2);
    % [M iM]=max(betaKernel,[],2);
%     iM=sortrows([iM (1:1:length(iM))'],1);
%     subplot(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))),i)
%     imagesc(regBeta(iM(:,2),:),[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    plot(regBeta);
    ylim([y_limits(1) y_limits(2)]);
    title(regLabels{betaIdx(i)});
    fig_configAxis(gca);
    offsetAxes(gca);
end
sgtitle(['Neuron #' num2str(id_neuron)])
