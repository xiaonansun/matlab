function dlc_plot_training_status(root_dir)

% %  projectDir is the root directory of the DLC project(s) 
% pDir=projectDir;
% csvFileName='learning_stats.csv';
% pDirList = split(genpath(pDir),';');
% csvFileList = cell(1,length(pDirList));
% 
% for i = 1:length(pDirList)
%     csvFileList(i) = dir(fullfile(pDirList{i},csvFileName));
% end

%% Input: provide a general root directory of where the training statistics file is located doesn't have to be precisely the directory of where it's located. Can even provide the drive, but this will take more time to locate the files.
% root_dir = '\\grid-hs\churchland_hpc_home\xisun\dlc_training'; % comment this line when executing as a function

file_name = 'learning_stats.csv';

filelist = dir(fullfile(root_dir, ['**\' file_name]));  %get list of files and folders in any subfolder

if ~isempty(filelist)
    [~,ind_sort] = sort(cell2mat({filelist.datenum}));
    path_most_recent = fullfile(filelist(ind_sort(end)).folder,filelist(ind_sort(end)).name);
    disp(['The most recent (' datestr(datetime(filelist(ind_sort(end)).datenum,'ConvertFrom','datenum')) ') DLC training analysis statistics will be plotted.'])
    training_stats = readtable(fullfile(filelist(ind_sort(end)).folder,filelist(ind_sort(end)).name));
else
    disp(['There are no DLC training learning statistics file (' file_name ') in any subdirectory'])
    return
end
str_dlc_model = 'dlc-models';
sub_dirs = regexp(path_most_recent,'\','split');
idx_dlc_model = find(contains(sub_dirs,str_dlc_model));

%% Plot
hLoss = figure; 
semilogy(training_stats{:,1},training_stats{:,2},'.k'); grid on;
% plot(training_stats{:,1},training_stats{:,3})
title([sub_dirs{idx_dlc_model-1} ',' sub_dirs{idx_dlc_model+1}])
ylabel('Loss (log_{10})'); xlabel('Iteration'); 
ax = gca;
fig_configAxis(ax);
fig_save_name = 'learning_stats_loss_vs_iteration';
fig_save_path = fullfile(filelist(ind_sort(end)).folder,fig_save_name);
saveas(hLoss,[fig_save_path '.png']);
saveas(hLoss,[fig_save_path '.eps'],'eps');
