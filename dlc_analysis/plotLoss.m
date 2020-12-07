function plotLoss(rootDir,compDir)

% Plots the latest DLC training status. rootDir is the project directory.
% compDir is the directory of the csv file for comparison; this directory
% can be a root directory

% dataFile = dir([rootDir filesep '**' filesep 'learning_stats.csv']);
dataFile = dir(fullfile(rootDir,'**','learning_stats.csv'));
[datetimeVal idx] = max(datenum({dataFile.date}));
table = readmatrix(fullfile(dataFile(idx).folder,dataFile(idx).name));
% table = table{:,:};

hFig1 = figure(1);

plot(table(:,1),log10(table(:,2)),'.-k');
title(['iter=' num2str(table(end,1)) '; Loss=' num2str(max(table(:,2))) ' at ' datestr(datetimeVal)]);
xlabel('Iteration');
ylabel('Loss (log_{10})');
% set(gca,'grid','on');

hold on;

if exist('compDir','var')
    dataFile = dir(fullfile(compDir,'**','learning_stats.csv'));
    [datetimeVal idx] = max(datenum({dataFile.date}));
    table = readmatrix(fullfile(dataFile(idx).folder,dataFile(idx).name));
    plot(table(:,1),log10(table(:,2)),'.-r');
end



%%

% rootDir = 'Y:\xisun\dlc_training\CSP27-Richard-2020-09-04';
% dataFile = dir(fullfile(rootDir,'**','*filtered.csv'));
% [datetimeVal idx] = max(datenum({dataFile.date}));
% colNames = readcell(fullfile(dataFile(idx).folder,dataFile(idx).name));
% pos=cell2mat(colNames(4:end,:));
% colNames(4:end,:)=[];
% 
% strTongue='Tongue'; strX='x'; strY='y'; strP='likelihood';
% idxColTongue=find(ismember(colNames(2,:),strTongue));
% x_tongue = idxColTongue(ismember(colNames(3,idxColTongue),strX));
% y_tongue =idxColTongue(ismember(colNames(3,idxColTongue),strY));
% p_tongue = idxColTongue(ismember(colNames(3,idxColTongue),strP));
% 
% 
% 
% hFig1 = figure(1);
% 
% plot(trajectories(:,1),log10(trajectories(:,2)),'.-k');
% title(['iter=' num2str(trajectories(end,1)) '; Loss=' num2str(max(trajectories(:,2))) ' at ' datestr(datetimeVal)]);
% xlabel('Iteration');
% ylabel('Loss (log_{10})');
% % set(gca,'grid','on');
