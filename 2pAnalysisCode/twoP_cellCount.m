rootdir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';

filelist = dir(fullfile(rootdir, '**\Vc.mat'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

for i = 1:length(filelist)
    Vc = load(fullfile(filelist(i).folder,filelist(i).name));
    Vc = Vc.Vc;
    redCell = readNPY(fullfile(filelist(i).folder,'redcell.npy'));
    isCell = readNPY(fullfile(filelist(i).folder,'iscell.npy'));
    sessionCounts{i,1} = size(Vc,1);
    sessionCounts{i,2} = size(Vc,3);
    sessionCounts{i,3} = sum(redCell(:,1).*isCell(:,1));
    sessionCounts{i,4} = filelist(i).folder;
end


%%
IndexC = strfind(sessionCounts(:,3),'CSP');
% cellTypeInd = cellfun(@(a)~isempty(a)&&a>0,IndexC);

Index = find(not(cellfun('isempty',IndexC)));

%%
cellCount = cell2mat(sessionCounts(:,1));
cellCountRed = cell2mat(sessionCounts(:,3));

mCellCount = mean(cellCount(Index));
mCellCountRed = mean(cellCountRed(Index));

stdCellCount = std(cellCount(Index));
stdCellCountRed = std(cellCountRed(Index));