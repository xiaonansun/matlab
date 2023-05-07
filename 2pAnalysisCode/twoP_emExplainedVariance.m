%% Compute cvR^2 of individual neurons
clear all;
animal = 'Plex50'; session = '200401a';

% lFile={'interpVc'; 'orgdimBeta';'orgfullcorr';'orgregData'};
% lFile={'interpVc';'interpVmotor';'interpVopMotor';'interpVtask';'orgVspontMotor';'orgfullcorr'};
lFile={'interpVc';'interpVspontMotor';'interpVopMotor';'interpVtask';'orgVspontMotor'};
%% 'orgfullcorr'
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
num_of_trials = size(Vc,2)/frames;
fmRsq = single(zeros(size(Vc,1),1));
uninsRsq = single(zeros(size(Vc,1),1));
insRsq = single(zeros(size(Vc,1),1));
taskRsq = single(zeros(size(Vc,1),1));

neu_tr_Rsq = single(zeros(size(Vc,1),num_of_trials));

for i = 1:size(Vc,1)
    full_cc = corrcoef(Vc(i,:),Vm(i,:)); 
    fmRsq(i) = full_cc(1,2)^2; % compute the R^2 for each neuron across all sessions% compute the correlation coefficient between the model versus data
    unins_cc = corrcoef(Vc(i,:),VspontMotor(i,:));
    uninsRsq(i) = unins_cc(1,2)^2; 
    ins_cc = corrcoef(Vc(i,:),VopMotor(i,:));
    insRsq(i) = ins_cc(1,2)^2;
    task_cc = corrcoef(Vc(i,:),Vtask(i,:));
    taskRsq(i) = task_cc(1,2)^2;
    
%     for j = 1:num_of_trials
%         trial_frames = frames*(j-1)+1:frames*j;
%         neu_tr_Vc = Vc(i,trial_frames);
%         neu_tr_Vm = Vm(i,trial_frames);
%         full_cc = corrcoef(neu_tr_Vc,neu_tr_Vm); % compute the single-neuron, single-trial correlation coefficient between the model versus data
%         neu_tr_Rsq(i,j) = full_cc(1,2)^2;
%     end
end

% neuron_rank = repmat((1:size(neu_tr_Rsq,1))',size(neu_tr_Rsq,2),1);
% trial_rank = repmat(1:size(neu_tr_Rsq,2),size(neu_tr_Rsq,1),1); trial_rank = trial_rank(:);
% [sorted_neu_tr_Rsq,idx_sorted_neu_tr_Rsq]=sort(neu_tr_Rsq(:),'descend');
% sorted_neuron_rank = neuron_rank(idx_sorted_neu_tr_Rsq);
% sorted_trial_rank = trial_rank(idx_sorted_neu_tr_Rsq);
% sorted_neu_tr_id = [sorted_neuron_rank sorted_trial_rank];
%%

hFig = figure(1); clf; hold on;
hFig.Position = [500 500 500 500];
x_data = 1:length(fmRsq);
[rankRsq,idxRsq] = sort(fmRsq,'descend');
line(x_data,rankRsq,'color',[0.5 0.5 0.5],'LineWidth',3);
line(x_data,uninsRsq(idxRsq),'color','k');
line(x_data,insRsq(idxRsq),'color','b');
line(x_data,taskRsq(idxRsq),'color','g');

xlabel('Neurons'); ylabel('cvR^2'); title(['Explained variance. ' animal ' ' session]);
hLeg = legend('Full Model','Uninstructed','Instructed','Task');
hLeg.Box = 'off'; hLeg.FontSize = 12;

fig_configAxis(gca);
offsetAxes(gca);

saveas(hFig,fullfile(dataDir,'exp_var_ind_neurons.pdf'));
saveas(hFig,fullfile(dataDir,'exp_var_ind_neurons.png'));
% figure
% imagesc(neu_tr_Rsq);
% colors = cbrewer('seq','Greys',1024);
% colormap(colors);
% xlabel('Trials'); ylabel('Neurons'); title('Single-neuron, single-trial cvR^2')