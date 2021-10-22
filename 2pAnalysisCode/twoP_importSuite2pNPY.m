function npy = twoP_importSuite2pNPY(animal, session)
% imports Fall.mat and chan2 data

% Requires s2ptomat.py
% Requires readNPY, part of npy-matlab package

% animal = 'Plex50'; % Input: Animal ID
% session = '200401b'; % Input: Session ID

% if ~exist('base_dir','var')
%     suite2p_base_dir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
% else 
%     suite2p_base_dir = base_dir;
% end

% suite2p_output_dir = fullfile(suite2p_base_dir, animal, 'imaging', session, 'suite2p\plane0');

% s2ptomat_path = '"C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab\2pAnalysisCode\s2ptomat.py"';
S = twoP_settings;

suite2p_output_dir = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir);
codeDir = fullfile(S.dir.codeRootDir,'2pAnalysisCode');
pyCodePath = fullfile(codeDir,'s2ptomat.py'); pyCodePath = ['"' pyCodePath '"'];
[status, result] = system(['python ' pyCodePath ' ' suite2p_output_dir]); % Converts ops.npy and stat.npy to ops.mat and stat.mat
load(fullfile(suite2p_output_dir, 'ops.mat'),'ops'); load(fullfile(suite2p_output_dir, 'stat.mat'),'stat');
npy.ops = ops{1,1}; npy.ops.animal = animal; npy.ops.session = session; npy.stat = stat; 
clear ops; clear stat;
npy_content = dir(fullfile([suite2p_output_dir], '*.npy')); npy_content = {npy_content.name};

for i_npy = 1:length(npy_content) % loads all npyfiles except for ops and stat, which cannot be read by readNPY.m
    [filepath, name, ext] = fileparts(fullfile(suite2p_output_dir, npy_content{i_npy}));
    if name ~= convertCharsToStrings('ops') && name ~= convertCharsToStrings('stat') && name ~= convertCharsToStrings('stat_orig')
        npy.(name) = readNPY(fullfile(suite2p_output_dir, npy_content{i_npy}));
    else
        continue
    end
end

iSub = strfind(suite2p_output_dir,filesep);
% subDirs = regexp(suite2p_output_dir,'\','split');
% [s2pPath,~,~] = fileparts(suite2p_output_dir);

bin_MScan_dir = suite2p_output_dir(1:iSub(end-1)); 
bin_MScan_dir_content = dir(bin_MScan_dir);
bin_MScan_dir_content = dir(fullfile(bin_MScan_dir,'*.bin'));

% bin_MScan_dir = fullfile(suite2p_base_dir, animal, 'imaging', session, filesep); 
% if size(bin_MScan_dir_content,1) > 1
    for i = 1:size(bin_MScan_dir_content,1)
        npy.bin_MScan_filepath{i} = fullfile(bin_MScan_dir_content(i).folder,bin_MScan_dir_content(i).name);
    end
% else
% npy.bin_MScan_filepath = [bin_MScan_dir bin_MScan_dir_content(contains({bin_MScan_dir_content.name}','bin')).name];
% end
% npy.bin_MScan_filepath = bin_MScan_filepath;

npy.bin_chan1_filepath = fullfile(suite2p_output_dir, 'data.bin');
npy.bin_chan2_filepath = fullfile(suite2p_output_dir, 'data_chan2.bin');

disp(['Number of cells: ' num2str(sum(npy.iscell(:,1))) ' | Total ROIs: ' num2str(numel(npy.iscell(:,1)))]);
disp(['Number of red cells: ' num2str(sum(npy.redcell(:,1)))]);
