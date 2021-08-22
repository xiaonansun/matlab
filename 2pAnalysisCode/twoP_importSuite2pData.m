function npy = twoP_importSuite2pData(animal, session, base_dir)
% dbstop
% imports Fall.mat and chan2 data

% Requires s2ptomat.py
% Requires readNPY, part of npy-matlab package
% All variable

% clear all;

% animal = 'Plex50'; % Input: Animal ID
% session = '200401b'; % Input: Session ID
if ~exist('base_dir','var')
    suite2p_base_dir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
else 
    suite2p_base_dir = base_dir;
end

suite2p_output_dir = [suite2p_base_dir filesep animal filesep 'imaging' filesep session filesep 'suite2p\plane0'];
s2ptomat_path = '"C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab\2pAnalysisCode\s2ptomat.py"';
[status, result] = system(['python ' s2ptomat_path ' ' suite2p_output_dir]); % Converts ops.npy and stat.npy to ops.mat and stat.mat
load([suite2p_output_dir filesep 'ops.mat']); load([suite2p_output_dir filesep 'stat.mat']);
npy.ops = ops{1,1}; npy.ops.animal = animal; npy.ops.session = session; clear ops;
npy.stat = stat; clear stat;
npy_content = dir(fullfile([suite2p_output_dir], '*.npy')); npy_content = {npy_content.name};

for i_npy = 1:length(npy_content) % loads all npyfiles except for ops and stat, which cannot be read by readNPY.m
    [filepath, name, ext] = fileparts([suite2p_output_dir filesep npy_content{i_npy}]);
    if name ~= convertCharsToStrings('ops') && name ~= convertCharsToStrings('stat')
        npy.(name) = readNPY([suite2p_output_dir filesep npy_content{i_npy}]);
    else
        continue
    end
end

bin_MScan_dir = [suite2p_base_dir filesep animal filesep 'imaging' filesep session filesep]; bin_MScan_dir_content = dir(bin_MScan_dir);
bin_MScan_filepath = [bin_MScan_dir bin_MScan_dir_content(contains({bin_MScan_dir_content.name}','bin')).name];
npy.bin_MScan_filepath = bin_MScan_filepath;

npy.bin_chan1_filepath = [suite2p_output_dir filesep 'data.bin'];
npy.bin_chan2_filepath = [suite2p_output_dir filesep 'data_chan2.bin'];

disp(['Number of cells: ' num2str(sum(npy.iscell(:,1))) ' | Total ROIs: ' num2str(numel(npy.iscell(:,1)))]);
disp(['Number of red cells: ' num2str(sum(npy.redcell(:,1)))]);
