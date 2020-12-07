session = '00038';
suite2p_base_dir = 'G:\2PData';
suite2p_output_dir = [suite2p_base_dir filesep session filesep 'suite2p\plane0'];
s2ptomat_path = '"C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab\2pAnalysis\s2ptomat.py"';
[status, result] = system(['python ' s2ptomat_path ' ' suite2p_output_dir]); % Converts ops.npy and stat.npy to ops.mat and stat.mat
load([suite2p_output_dir filesep 'ops.mat']); load([suite2p_output_dir filesep 'stat.mat']);
npy.ops = ops{1,1}; npy.ops.session = session; clear ops;
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

%%

cellIdx=find(npy.iscell(:,1));
plot(npy.spks(cellIdx(2),:),npy.spks(cellIdx(212),:),'.k');
r = xcorr(npy.spks(cellIdx(150),:),npy.spks(cellIdx(212),:));
pearson_corr = corr([npy.spks(cellIdx(150),:),npy.spks(cellIdx(212),:)],'Type','spearman');
[r,p]=corrcoef(npy.spks(cellIdx(151),:),npy.spks(cellIdx(212),:));

ccMat=zeros(length(cellIdx));
ccMatSorted=zeros(length(cellIdx));
ccMatSortedIdx=zeros(length(cellIdx));

% ccMatSorted=ccMat;
tic
parfor i = 1:length(cellIdx)
    v = zeros(1, length(cellIdx));
    for j = 1:length(cellIdx)
       [r,p]=corrcoef(npy.spks(cellIdx(i),:),npy.spks(cellIdx(j),:));
       v(j)=r(1,2);
    end
    ccMat(i, :) = v;
%     [B,I]=sort(ccMat(i,1+i:end),2);
%     ccMatSorted(i,:)=I;
%     ccMatSorted(i+1:end,i)=I;
end
toc


for i = 1:length(cellIdx)
    [B,I]=sort(ccMat(i,1+i:end),2);
    ccMatSortedIdx(i,1+i:end)=fliplr(I);
    t=ccMat(i,1+i:end);
    ccMatSorted(i,1+i:end) = t(ccMatSortedIdx(i,1+i:end));
    
%     ccMatSorted(i+1:end,i)=I;
end
ccMatSorted= ccMatSorted+fliplr(rot90(ccMatSorted,-1));

imagesc(ccMatSorted);
% [B,I]=sort(ccMat(1,2:end),2);

% imagesc(ccMat);

% [B,I]=sort(ccMat,1);

% imagesc(ccMat(I(:,1),I(1,:)))

% 
% bin_MScan_dir = [suite2p_base_dir filesep animal filesep session filesep]; bin_MScan_dir_content = dir(bin_MScan_dir);
% bin_MScan_filepath = [bin_MScan_dir bin_MScan_dir_content(contains({bin_MScan_dir_content.name}','bin')).name];
% npy.bin_MScan_filepath =bin_MScan_filepath;
% 
% npy.bin_chan1_filepath = [suite2p_output_dir filesep 'data.bin'];
% npy.bin_chan2_filepath = [suite2p_output_dir filesep 'data_chan2.bin'];
% 
% disp(['Number of cells: ' num2str(sum(npy.iscell(:,1))) ' | Total ROIs: ' num2str(numel(npy.iscell(:,1)))]);
% disp(['Number of red cells: ' num2str(sum(npy.redcell(:,1)))]);