%%
clear all;
data_path = "H:\lukas\miniscope_data.mat";
[data_dir,~,~] = fileparts(data_path);
% [ridgeVals, dimBeta, Vm, fullBeta, fullR, fullIdx, fullRidge, fullLabels, fullMap, fullMovie] = choiceModel_small(data_path);

% function [ridgeVals, dimBeta, Vm, fullBeta, fullR, fullIdx, fullRidge, fullLabels, fullMap, fullMovie] = choiceModel_small(iDataPath)
% Load data
load(data_path,'frames_per_trial','reg_group_name','reg_idx','X','Y');
%%

% Vc is the neural data
% fullR is the design matrix
%% Cross validation
% [Vm, fullBeta, fullR, fullIdx, fullRidge, fullLabels, fullMap, fullMovie] = crossValModel(regLabels); % full means both task and motor regressors are included in the cross validation
% regLabels = cellstr(reg_group_name);
f = @(a) regexprep(a,' ','_'); % define cellfun to replace space with underscore
regLabels = cellfun(f,cellstr(reg_group_name),'UniformOutput',false); % replace space to underscore in regLabels
regIdx = reg_idx+1; % python indexing starts at 0, while matlab starts at 1. reg_idx is from python data
fullR = X; Vc = Y'; % Initializes some variables as input to subsequent functions
motorLabels = {'Time' 'Head_direction' 'DLC' 'Video_and_ME_components'};
taskLabels = {'Time' 'Cognitive'};
opMotorLabels = {'Time' 'Head_direction' 'DLC'};
spontMotorLabels  = {'Time' 'Video_and_ME_components'};
frames = frames_per_trial;

%% Full and motor model model ridge regression
[ridgeVals, dimBeta] = ridgeMML(Y, X, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for full model: %f\n', mean(ridgeVals));
fullVm = (fullR * dimBeta)';
save(fullfile(data_dir,'orgdimBeta.mat'), 'dimBeta', 'ridgeVals');

mInd = ismember(regIdx(~rejIdx), find(ismember(regLabels,motorLabels)));
motorR = fullR(:, mInd);
[motorRidge, motorBeta] = ridgeMML(Y, motorR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original dataset, motor only model: %f\n', mean(motorRidge));
Vm = (motorR * motorBeta)';
save(fullfile(data_dir,'orgVMotor.mat'), 'Vm', 'frames'); %save predicted data based on motor model

mInd = ismember(regIdx(~rejIdx), find(ismember(regLabels,motorLabels(~ismember(motorLabels,opMotorLabels)))));
spontMotorR = fullR(:, mInd);
[spontMotorRidge, spontMotorBeta] = ridgeMML(Y, spontMotorR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original dataset, spont-motor only model: %f\n', mean(spontMotorRidge));
Vm = (spontMotorR * spontMotorBeta)';
save(fullfile(data_dir,'orgVspontMotor.mat'), 'Vm', 'frames'); %save predicted data based on spont motor model

%% Reduced model ridge regression
num_reduced = 1;
clear fullCorr
fullCorr = zeros(1,size(Vc,1));
for x = 1 : size(Vc,1)
    fullCorr(x) = corr2(Vc(x,:),fullVm(x,:)); % get correlation between data and model
end
dRsq = cell(length(regLabels)-num_reduced,1);
redCorr = cell(length(regLabels)-num_reduced,1);
for i = (1+num_reduced):length(regLabels) % cycles through the regressors to be removed for the reduced model (time is never removed)
    %%
    redLabels = regLabels;
    redLabels(i) = []; % removes a one regressor at a time
    redInd = ismember(regIdx(~rejIdx), find(ismember(regLabels,redLabels))); % generates the indices of the reduced model
    redR = fullR(:, redInd); % Design matrix where one regressor is removed
    [redRidge, redBeta] = ridgeMML(Y, redR, true); %get ridge penalties and beta weights.
%     fprintf(['Mean ridge penalty for ' regLabels{i} ' model: %f\n', mean(motorRidge)]);
    Vm = (redR * redBeta)';
    for x = 1 : size(Vc,1)
        redCorr{i-num_reduced}(x) = corr2(Vc(x,:),Vm(x,:)); % get correlation between data and reduced model
    end
    dRsq{i-1} = fullCorr.^2-redCorr{i-num_reduced}.^2; % computes deltaR^2
end

%% Plot reduced model
hFig_dRsq = figure(2); clf; hold on;
hFig_dRsq.Position = [500 500 500 500];
cell_id = 1:size(Y,2);
[~,sorted_idx] = sort(fullCorr,'descend');
CT=cbrewer('qual', 'Dark2', 8);
for i = 1:length(dRsq)
line(cell_id,dRsq{i}(sorted_idx).^2,'Color',CT(i,:),...
    'LineWidth',2);
end

xlabel('Neurons'); ylabel('\DeltaR^2'); title(['Explained variance']);
hLeg = legend(regLabels(2:end),'Interpreter','none');
hLeg.Box = 'off'; hLeg.FontSize = 12;

fig_configAxis(gca);
offsetAxes(gca);

saveas(hFig_dRsq,fullfile(data_dir,'dRsq_ind_neurons.pdf'));
saveas(hFig_dRsq,fullfile(data_dir,'dRsq_ind_neurons.png'));

%% Cross validation
% [Vm, cBeta, cR, subIdx, cRidge, cLabels, cMap] =  crossValModel_unnested(cLabels,regLabels,regIdx,rejIdx,Vc,fullR)

[Vm, fullBeta, fullR, fullIdx, fullRidge, fullLabels, fullMap] =  crossValModel_unnested(regLabels,regLabels,regIdx,rejIdx,Y',X);
save(fullfile(data_dir,'orgfullcorr.mat'), 'Vm', 'fullBeta', 'fullIdx', 'fullR', 'fullLabels', 'fullRidge', 'regLabels', 'fullMap','-v7.3');

[Vmotor, motorBeta, motorR, motorIdx, motorRidge, motorLabels, motorMap] = crossValModel_unnested(motorLabels,regLabels,regIdx,rejIdx,Y',X);
fprintf('Mean ridge penalty for motor-only, zero-mean model: %f\n', mean(motorRidge));
save([fullfile(data_dir,'interpVmotor.mat')], 'Vmotor', 'frames'); %save predicted data based on motor model
save([fullfile(data_dir,'motorBeta.mat')], 'motorBeta', 'motorRidge');
save(fullfile(data_dir,'motorregData.mat'), 'motorR', 'motorIdx', 'motorLabels', 'motorMap','-v7.3');

[Vtask, taskBeta, taskR, taskIdx, taskRidge, taskLabels, taskMap] = crossValModel_unnested(taskLabels,regLabels,regIdx,rejIdx,Y',X);
fprintf('Mean ridge penalty for task-only, zero-mean model: %f\n', mean(taskRidge));
save(fullfile(data_dir,'interpVtask.mat'), 'Vtask', 'frames'); %save predicted data based on task model
save(fullfile(data_dir,'taskBeta.mat'), 'taskBeta', 'taskRidge');
save(fullfile(data_dir,'taskregData.mat'), 'taskR', 'taskIdx', 'taskLabels', 'taskMap', '-v7.3');

[VopMotor, opMotorBeta, opMotorR, opMotorIdx, opMotorRidge, opMotorLabels, opMotorMap] = crossValModel_unnested(opMotorLabels,regLabels,regIdx,rejIdx,Y',X);
fprintf('Mean ridge penalty for opMotor-only, zero-mean model: %f\n', mean(opMotorRidge));
save(fullfile(data_dir,'interpVopMotor.mat'), 'VopMotor', 'frames'); %save predicted data based on operant motor model
save(fullfile(data_dir,'opMotorBeta.mat'), 'opMotorBeta', 'opMotorRidge');
save(fullfile(data_dir,'opMotorregData.mat'), 'opMotorR', 'opMotorIdx', 'opMotorLabels', 'opMotorMap', '-v7.3');


[VspontMotor, spontMotorBeta, spontMotorR, spontMotorIdx, spontMotorRidge, spontMotorLabels, spontMotorMap] = crossValModel_unnested(spontMotorLabels,regLabels,regIdx,rejIdx,Y',X);
fprintf('Mean ridge penalty for spontMotor-only, zero-mean model: %f\n', mean(spontMotorRidge));
save(fullfile(data_dir,'interpVspontMotor.mat'), 'VspontMotor', 'frames'); %save predicted data based on spontaneous motor model
save(fullfile(data_dir,'spontMotorBeta.mat'), 'spontMotorBeta', 'spontMotorRidge');
save(fullfile(data_dir,'spontMotorregData.mat'), 'spontMotorR', 'spontMotorIdx', 'spontMotorLabels', 'spontMotorMap', '-v7.3');

%% Plots the cross-validated R^2 
hFig_Rsq = figure(2); clf; hold on;
hFig_Rsq.Position = [500 500 500 500];
cell_id = 1:size(Y,2);
[~,sorted_idx] = sort(fullMap,'descend');
line(cell_id,fullMap(sorted_idx).^2,'Color',[0.5 0.5 0.5],...
    'LineWidth',4);
line(cell_id,taskMap(sorted_idx).^2,'Color','g');
line(cell_id,opMotorMap(sorted_idx).^2,'color','b');
line(cell_id,spontMotorMap(sorted_idx).^2,'color','k');

xlabel('Neurons'); ylabel('cvR^2'); title(['Explained variance']);
hLeg = legend('Full Model','Task','Instructed','Uninstructed');
hLeg.Box = 'off'; hLeg.FontSize = 12;

fig_configAxis(gca);
offsetAxes(gca);

saveas(hFig_Rsq,fullfile(data_dir,'exp_var_ind_neurons.pdf'));
saveas(hFig_Rsq,fullfile(data_dir,'exp_var_ind_neurons.png'));