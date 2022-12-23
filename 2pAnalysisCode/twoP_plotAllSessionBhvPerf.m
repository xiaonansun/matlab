function twoP_plotAllSessionBhvPerf(animal)

%%
% animal='CSP30';
min_file_size = 1000000; % minimum file size
S = twoP_settings;
file_dir = fullfile(S.dir.bhvRootDir,animal,S.dir.bhvSubDir);
save_analysis_dir = fullfile(S.dir.imagingRootDir,'behavior_analysis');
if ~exist(save_analysis_dir,'dir'); mkdir(save_analysis_dir); end
file_list = dir(file_dir);
file_list = file_list(cell2mat({file_list.bytes})> min_file_size);
for i = 1:length(file_list) % Load the first-trial start time of each session from the behavior file
    echo off
    load(fullfile(file_list(i).folder,file_list(i).name),'SessionData');
    dt = twoP_displaySessionTime(SessionData);
    trial_start_time(i) = dt{1};
end

% file_name_str = {file_list_by_date.name};
% file_date_str = cellfun(@(x) x,trial_start_time,'UniformOutput',false);
% file_date_str = datetime(trial_start_time,'InputFormat','MM/dd/yy HH:mm:ss.SSS');
[~,idx_date] = sort(trial_start_time);
% [~,idx_date] = sort(cell2mat({file_list.datenum}));

% for i = 1:length(idx_date)
%     file_list_by_date(i) = file_list(idx_date(i));
% end

file_list_by_date = file_list(idx_date);
%%
hFig = figure;
set(gcf,'Units','inches',...
    'Position',[0 0 8.5 8.5]);
hTile = tiledlayout(ceil(sqrt(length(file_list))),ceil(sqrt(length(file_list))),...
    'TileSpacing','tight',...
    'Padding','tight');
for i = 1:length(file_list_by_date)
    nexttile;
    load(fullfile(file_list_by_date(i).folder,file_list_by_date(i).name),'SessionData');
    [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = rateDisc_audioDiscCurve(SessionData,1:length(SessionData.Rewarded));
    plot(distRatio,pChoseHigh,'.-k');
    ax = fig_configAxis(gca);
    fig_separateAxes;
    xlim([0 1]); ylim([0 1]);
    hTitle = title([file_list_by_date(i).name(end-22:end-13) '_' file_list_by_date(i).name(end-4)],...
        'Interpreter','none',...
        'FontSize',6);
    if i == 1 || rem(i,ceil(sqrt(length(file_list)))) == 1
        ylabel('p_{right}');
    else
        set(ax,'YTickLabel',[]);
    end
    if i > length(file_list)-ceil(sqrt(length(file_list)))
        xlabel('F_{right}');
    else
        set(ax,'XTickLabel',[]);
    end
    
end

title(hTile,animal)
save_file_path = fullfile(save_analysis_dir,[animal '_single_session_performance.pdf']);
saveas(hFig,save_file_path);
