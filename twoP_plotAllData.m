% Before running this script, twoP_combineData.m should have been ran
exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1); colExpertise = exps(:,5); colSession = exps(:,6);

S = twoP_settings;

if S.isUnix == true; parpool('local',32); end

celltype = cell(size(exps,1),1);
sCellType = contains(colAnimal,S.cellTypes);
sExpertise = contains(colExpertise,S.expertise(2:end));
s = find(sCellType & sExpertise);
c = zeros(length(s),1);
D = cell(size(exps,1),1);
E = cell(size(exps,1),1);
CT = cell(size(exps,1),1);

parfor j = 1:length(s)
    try
        animal = colAnimal{s(j)}; session = colSession{s(j)};
        npyPathIsCell = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'iscell.npy');
        npyPathRedCell = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'redcell.npy');
        idxIC = readNPY(npyPathIsCell); idxRC = readNPY(npyPathRedCell);
        idxRC = idxRC(logical(idxIC(:,1)),1);
        pathVc = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'Vc.mat');
        temp = load(pathVc,'Vc'); Vc = temp.Vc;  temp = [];
        pathcBhv = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'cBhv.mat');
        temp = load(pathcBhv,'cBhv'); cBhv = temp.cBhv;  temp=[];
        for i = 1:length(S.cellTypes)
            if contains(colAnimal{s(j)},S.cellTypes{i})
                CT{j} = S.cellTypes{i};
                break
            else
                continue
            end
        end
        D{j} = mean(Vc,3,'omitnan');
        E{j} = idxRC;
        c(j) = 1;
        disp(['Loaded ' animal '_' session]);
    end
end
%%

msPerFrame=32.3638; sPerFrame = msPerFrame/1000;
frameRate = 1000/msPerFrame; % Frame rate of imaging
sRate = frameRate; segFrames = cumsum(floor(S.segIdx * sRate)); % max nr of frames per segment
sF = [1 segFrames(2:end)];
idxOnset = sF(2);
% xTL = (-sF(2)+1)*sPerFrame:sPerFrame:(sF(end)-sF(2))*sPerFrame;
lineWidth = 2;

idxD = ~cellfun(@isempty,D);
A = vertcat(D{idxD}); irA = vertcat(E{idxD}); % ir = index of red cells;
NR = A(~irA,:);
CSP = vertcat(D{idxD & strcmp(CT,'CSP')}); irCSP = vertcat(E{idxD & strcmp(CT,'CSP')});
Plex = vertcat(D{idxD & strcmp(CT,'Plex')}); irPlex = vertcat(E{idxD & strcmp(CT,'Plex')});
Fez = vertcat(D{idxD & strcmp(CT,'Fez')}); irFez = vertcat(E{idxD & strcmp(CT,'Fez')});

mA = mean(A,1,'omitnan');
mNR = mean(NR,1,'omitnan');
mCSP = mean(CSP(logical(irCSP),:),1,'omitnan');
mPlex = mean(Plex(logical(irPlex),:),1,'omitnan');
mFez = mean(Fez(logical(irFez),:),1,'omitnan');
x = 1:size(mA,2); t = (-sF(2)+1)*sPerFrame:sPerFrame:(sF(end)-sF(2))*sPerFrame;

fPETH = figure(1);
set(fPETH, 'Position',[250 250 500 200]);
lAll = line(t,twoP_normToOnset(mA),'Color','r',...
    'LineWidth',lineWidth);
lNR = line(t,twoP_normToOnset(mNR),'Color','k',...
    'LineWidth',lineWidth);
lCSP = line(t,twoP_normToOnset(mCSP),'Color',[0.9290 0.6940 0.1250],...
    'LineWidth',lineWidth);
lPlex = line(t,twoP_normToOnset(mPlex),'Color','g',...
    'LineWidth',lineWidth);
lFez = line(t,twoP_normToOnset(mFez),'Color','b',...
    'LineWidth',lineWidth);
xlineLabel = {'Handles','Stimulus','Delay','Response',' '};

xline(t(sF),'-k',xlineLabel,...
    'HandleVisibility','off',...
    'LabelVerticalAlignment','bottom');
legend({'All PyNs';'No tdT';'CSP';'PlexinD1';'Fezf2'},...
    'Color','none',...
    'Box','off',...
    'Location','Northwest',...
    'FontSize',6,...
    'NumColumns',2);

ax = gca;
ax = fig_configAxis(ax);

xlabel('Time from stimulus onset(s)'); ylabel({'Inferred spikes';'(normalized)'});
title('PETH - all sessions and all trials - by cell type');
exportgraphics(fPETH,fullfile(S.dir.imagingRootDir,'analysis','PETH_by_cell type.pdf'));


%%
A = cell(size(exps,1),1);
A{s(find(c))} = D{find(c)};