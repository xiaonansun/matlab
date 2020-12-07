%% imports Fall.mat and chan2 data
clear all;

animal = 'Plex50'; % Input: Animal ID
session = '200401b'; % Input: Session ID

suite2p_base_dir = 'W:\DOMAIN=CSHL\smusall';
suite2p_output_dir = [suite2p_base_dir filesep animal filesep session filesep 'suite2p\plane0'];
bin_dir = [suite2p_base_dir filesep animal filesep session filesep]; bin_dir_content = dir(bin_dir);
bin_filepath = [bin_dir bin_dir_content(contains({bin_dir_content.name}','bin')).name];
suite2p_Fall_path = [suite2p_output_dir filesep 'Fall.mat'];
load(suite2p_Fall_path);

binary_chanG_file_path = [suite2p_output_dir filesep 'data.bin'];
binary_chanR_file_path = [suite2p_output_dir filesep 'data_chanR.bin'];

Fneu_chan2 = readNPY([suite2p_output_dir filesep 'Fneu_chan2.npy']);
F_chan2 = readNPY([suite2p_output_dir filesep 'F_chan2.npy']);

disp(['Number of cells: ' num2str(sum(iscell(:,1))) ' | Total ROIs: ' num2str(numel(iscell(:,1)))]);
disp(['Number of red cells: ' num2str(sum(redcell(:,1)))]);
cell_idx = find(iscell(:,1)); isnotcell_idx = find(~iscell(:,1)); 
redcell_idx = find(redcell(:,1)); isnotredcell_idx = find(~redcell(:,1));


%% 
clf;

idx_ROI = 90;
num_of_bins = 400;
nCols = 3;
nRows = 4;
nPlots = nCols*nRows;
fSession = 0.05; %fraction of session to smooth
fps = 30.9;
time = 1/fps:1/fps:length(F)/fps;
nFramesSmooth = length(F)*fSession;
smoothMethod = 'movmean';
F_smooth = smoothdata(F,2,smoothMethod,nFramesSmooth);
Fneu_smooth = smoothdata(Fneu,2,smoothMethod,nFramesSmooth);
dF = F-F_smooth;
dFneu = Fneu-Fneu_smooth;
pdF_ROI = fitdist(dF(cell_idx(idx_ROI),:)','Normal');
pdFneu_ROI = fitdist(dFneu(cell_idx(idx_ROI),:)','Normal');
smoothLineWidth = 2;
smoothLineColor = 'r';

hFig1 = figure(1);
subplot(nRows,nPlots/nRows,1); hold on; 
title('F and F_{neu}');
histogram(F(cell_idx(idx_ROI),:),num_of_bins,'FaceColor','k'); 
histogram(Fneu(cell_idx(idx_ROI),:),num_of_bins, 'FaceColor','g','EdgeColor','none');
xlim([0 10000]);

subplot(nRows,nPlots/nRows,[2 3]); hold on;
title('Raw traces');
hF = plot(time,F(cell_idx(idx_ROI),:),'b','LineWidth',1); hF.Color(4) = 0.1;
hF_smooth = plot(time, F_smooth(cell_idx(idx_ROI),:),'b','LineWidth',2); 
hFneu = plot(time, Fneu(cell_idx(idx_ROI),:),'g','LineWidth',1); hFneu.Color(4) = 0.2;
hFneu_smooth = plot(time, Fneu_smooth(cell_idx(idx_ROI),:),'g','LineWidth',2);  hFneu.Color(4) = 0.2;
xlim([0 max(time)]);

subplot(nRows,nPlots/nRows,4); hold on;
title('Normalized F and F_{neu}');
hHist3 = histogram(dF(cell_idx(idx_ROI),:),num_of_bins);
hHist4 = histogram(dFneu(cell_idx(idx_ROI),:),num_of_bins,'FaceColor','g','EdgeColor','none','FaceAlpha',1);
% xlim([0 10000]);
subplot(nRows,nPlots/nRows,[5 6]); hold on; 
title('Normalized traces');
plot(time, dF(cell_idx(idx_ROI),:)); 
plot(time, dFneu(cell_idx(idx_ROI),:),'g'); 
xlim([0 max(time)]);

subplot(nRows,nPlots/nRows,7:9); hold on;
title('Standard deviation trace'); ylabel('sigma')
plot(time,dF(cell_idx(idx_ROI),:)/pdF_ROI.sigma); 
hFneu_sigma = plot(time, dFneu(cell_idx(idx_ROI),:)/pdF_ROI.sigma,'g'); hFneu_sigma.Color(4) = 0.2;
xlim([0 max(time)]);

subplot(nRows,nPlots/nRows,10); histogram(spks(cell_idx(idx_ROI),find(spks(cell_idx(idx_ROI),:))),num_of_bins);
subplot(nRows,nPlots/nRows,[11 12]); plot(time,spks(cell_idx(idx_ROI),:)); 
xlim([0 max(time)]);


% hFig2 = figure(2); hold on;
% plot(time,spks(cell_idx(idx_ROI),:)); 
% plot(time,medfilt1(spks(cell_idx(idx_ROI),:)),'g'); 

% a = find(spks(cell_idx(idx_ROI),:));

% x_values = floor(min(dF(cell_idx(idx_ROI),:))):1:ceil(max(dF(cell_idx(idx_ROI),:))); 
% F_pdf= pdf(pdF_ROI,x_values);
% hFig2 = figure(2); hold on;
% histogram(dF(cell_idx(idx_ROI),:),num_of_bins,'Normalization','probability');
% plot(x_values,F_pdf,'r','LineWidth',3);
%% import binary and imaging data using Simon's script

% data = twoP_alignDetectionTask(suite2p_Fall_path, bin_filepath);
[vals, sampFreq] = readMOMAnalog(bin_filepath);

%% Plot MScan binary data
plot_duration = 50; % duration to plot in seconds

figure(2); hold on;
for i = 1:size(vals,1)
    subplot(size(vals,1),1,i); xticks([]);
    plot(vals(i,1:plot_duration*sampFreq));
end

%%
data_samp = vals(2,1:plot_duration*sampFreq);
[pks1 frameOnset_idx1] = findpeaks(vals(2,:));
[pks2 frameOnset_idx2] = findpeaks(-vals(2,:));
frameOnset_idx = sort([frameOnset_idx1 frameOnset_idx2]); %% this is the index of the binary data when the mirror changes direction, which is the onset of a new frame
disp(['Total number of frames: ' num2str(length(frameOnset_idx))]);


%%
n_frames = 100;

Fneu_begin_end = [median(Fneu(cell_idx,1:n_frames),2) median(Fneu(cell_idx,n_frames:end),2)];
Fneu_delta = (Fneu_begin_end(:,2)-Fneu_begin_end(:,1))./Fneu_begin_end(:,1);

F_begin_end = [median(F(cell_idx,1:n_frames),2) median(F(cell_idx,n_frames:end),2)];
F_delta = (F_begin_end(:,2)-F_begin_end(:,1))./F_begin_end(:,1);

subplot(3,2,1); histogram(Fneu_delta,50);
subplot(3,2,3); histogram(F_delta,50);

subplot(3,2,2), plot(Fneu_delta,F_delta,'.k');

plot([1 2],Fneu_begin_end);

%% 
spks_max = max(spks(cell_idx,:),[],2);
hist(spks_max,50)
% set(gca,'YScale','log');

scatter(iscell(:,1),iscell(:,2))

%% load NPY
% npy_dir = [suite2p_file_dir filesep 'suite2p\plane0'];
% npy_file_struct = dir(fullfile(npy_dir,'*.npy'));
% npy_file_list = {npy_file_struct.name};
% 
% for i = 1:numel(npy_file_list)
%     readNPY([npy_dir filesep npy_file_list{i}]);
% end
%% load binaries
outputFPS = 120;
avg_frames = 30;
idx_first_frame = 10000;
num_of_frames_to_load = 1500;

bytes_per_pixel = 2; x_res = 512; y_res = 512;
bytes_per_frame = x_res*y_res*bytes_per_pixel;
bytes_per_line = x_res*bytes_per_pixel;

bin_frame_idx = bytes_per_frame*(idx_first_frame-1);

chanG_fileID = fopen(bin_filepath, 'r');
chanR_fileID = fopen(bin_chan2_filepath, 'r');

fseek(chanG_fileID,(idx_first_frame-1)*bytes_per_frame,'bof');
chanG = fread(chanG_fileID, [1 x_res*y_res*num_of_frames_to_load],'uint16');
chanG = reshape(chanG,x_res,y_res,[], size(chanG,2)/(x_res*y_res));
chanG_mean = movmean(chanG,avg_frames,4);

fseek(chanR_fileID,(idx_first_frame-1)*bytes_per_frame,'bof');
chanR = fread(chanR_fileID, [1 x_res*y_res*num_of_frames_to_load],'uint16');
chanR = reshape(chanR,x_res,y_res,[], size(chanR,2)/(x_res*y_res));
chanR_mean = movmean(chanR,avg_frames,4);

vid_output_G = rescale(chanG_mean); color_range = [0.3 0.8];
vid_output_G = rescale(vid_output_G,color_range(1),color_range(2));
vid_output_R = rescale(chanR_mean); color_range = [0.3 0.8];
vid_output_R = rescale(vid_output_R,color_range(1),color_range(2));

vid_output_both_chan = [vid_output_G vid_output_R];

saveFileName = datestr(datetime,'yyyy-mm-dd_HHMMSS');
[saveDir, fileNameG, ext] = fileparts(bin_filepath);
[saveDir, fileNameR, ext] = fileparts(bin_chan2_filepath);

vidObj = VideoWriter([saveDir filesep fileNameG '_' saveFileName '.avi']);
vidObj.FrameRate = outputFPS;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj, vid_output_G);
close(vidObj);

vidObjR = VideoWriter([saveDir filesep fileNameR '_' saveFileName '.avi']);
vidObjR.FrameRate = outputFPS;
vidObjR.Quality = 100;
open(vidObjR);
writeVideo(vidObjR, vid_output_R);
close(vidObjR);

vidObjGR = VideoWriter([saveDir filesep '_GR_' saveFileName '.avi']);
vidObjGR.FrameRate = outputFPS;
vidObjGR.Quality = 100;
open(vidObjGR);
writeVideo(vidObjGR, vid_output_both_chan);
close(vidObjGR);

handle = implay(vid_output_G,outputFPS);
handle.Visual.ColorMap.UserRange = 1;
handle.Visual.ColorMap.UserRangeMin = min(vid_output_G(:)-0.2);
handle.Visual.ColorMap.UserRangeMax = max(vid_output_G(:)+0.2);
%% Generate a mean image

frewind(chanG_fileID); idx = 1; chanG_sum = zeros(1,x_res*y_res);
while ~feof(chanG_fileID) && idx <= num_of_frames_to_load
chanG_sum = chanG_sum + fread(chanG_fileID, [1 x_res*y_res],'uint16');
fprintf(num2str(idx)); pause(0.01); clc; 
idx = idx+1;
end
chanG_mean = reshape(chanG_sum,x_res,y_res,[], size(chanG_sum,2)/(x_res*y_res))./idx;

%% testing fopen fread, and fseek 
fseek(chanG_fileID, 13,'bof');
F_pixel = fread(chanG_fileID,1,'uint16')

%% plot basic statistics
% hIsCell = histogram(iscell(:,2),50);
% hRedCell = histogram(redcell(:,2),100);
% set(gca,'YScale','log');
% set(hCell, 'FaceColor','k');

is_red_cell_idx = find(redcell(:,1));
disp(['Red cell probability: ' num2str(mean(redcell(is_red_cell_idx,2))) ' (mean); ' num2str(median(redcell(is_red_cell_idx,2))) ' (median); [' num2str(min(redcell(is_red_cell_idx,2))) ' ' num2str(max(redcell(is_red_cell_idx,2))) '] [min max]']);
% histogram(redcell(is_red_cell_idx,2),40);

not_red_cell_idx = find(~redcell(:,1));
disp(['Not red cell probability: ' num2str(mean(redcell(not_red_cell_idx,2))) ' (mean); ' num2str(median(redcell(not_red_cell_idx,2))) ' (median); [' num2str(min(redcell(not_red_cell_idx,2))) ' ' num2str(max(redcell(not_red_cell_idx,2))) '] [min max]']);

% histogram(redcell(not_red_cell_idx,2),40);