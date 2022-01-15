function [combData,meta]= twoP_loadSessions(cellType,expertise,location,depth)
%(1) cellType: CSP vs Plex vs Fez
%(2) expertise: Naive vs Trained
%(3) location: ALM vs MM
%(4) depth: Superficial vs Intermediate vs Deep

%%
% clear cellType expertise location

S = twoP_settings;
imagingRootDir = S.dir.imagingRootDir; imagingSubDir = S.dir.imagingSubDir;

if ~exist('cellType','var') || isempty(cellType)
cellType = {'CSP','Plex','Fez'};
end

if ~exist('expertise','var') || isempty(expertise)
expertise = {'Trained'};
end

if ~exist('location','var') || isempty(location)
location = {'ALM'};
end

if ~exist('depth','var') || isempty(depth)
depth = {'Intermediate'};
end

exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1); 
colDates = exps(:,2);
colLocation = exps(:,3); 
colDepth = exps(:,4);
colDepthCat = twoP_getImagingDepth(colDepth);
colExpertise = exps(:,5); 
colSession = exps(:,6);

idxCellType = contains(colAnimal,cellType);
if length(location) > 1 && iscell(location)
    idxLocation = contains(colLocation,location);
else
    idxLocation = strcmp(colLocation,location);
end

if length(expertise) > 1 && iscell(depth)
    idxExpertise = contains(colExpertise,expertise);
else
    idxExpertise = strcmp(colExpertise,expertise);
end

if length(depth) > 1 && iscell(depth)
    idxDepth = contains(colDepthCat,depth);
else
    idxDepth = strcmp(colDepthCat,depth);
end

idxExps = idxCellType & idxLocation & idxExpertise & idxDepth;

allAnimal = cell(size(exps,1),1);
allSession = cell(size(exps,1),1);
allDates = cell(size(exps,1),1);
allVc = cell(size(exps,1),1);
allVcNan = cell(size(exps,1),1);
allVcNormSD = cell(size(exps,1),1); % Inferred spiking activity normalized to units of standard deviation
allcBhv = cell(size(exps,1),1);
idxRedCell = cell(size(exps,1),1);

parfor i = 1:length(idxExps)
    %%
    try
        if idxExps(i) == 1
        Vc = load(fullfile(imagingRootDir,colAnimal{i},'imaging',colSession{i},imagingSubDir,'Vc.mat'),'Vc'); Vc = Vc.Vc;
        cBhv = load(fullfile(imagingRootDir,colAnimal{i},'imaging',colSession{i},imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
        idxCell = readNPY(fullfile(imagingRootDir,colAnimal{i},'imaging',colSession{i},imagingSubDir,'iscell.npy'));
        idxRed = readNPY(fullfile(imagingRootDir,colAnimal{i},'imaging',colSession{i},imagingSubDir,'redcell.npy'));
        idxRedCell{i} = idxRed(find(idxCell(:,1)));
        allAnimal{i} = colAnimal{i};
        allDates{i} = colDates{i};
        allSession{i} = colSession{i};
        allVc{i} = Vc;
        allVcNan{i} = sum(isnan(Vc),3);
        allVcNormSD{i} = bsxfun(@rdivide,Vc,std(Vc,0,3,'omitnan'));
        allcBhv{i} = twoP_bhvSubSelection(cBhv);
        disp(['Loaded ' allAnimal{i} ' ' allSession{i} '.']);
        end
    catch ME
    end
end

combData = [allAnimal allSession allVc allVcNormSD allcBhv idxRedCell allVcNan allDates];
meta.colLabels = {'Animal','Session','Vc','Z-score','cBhv','idxRedCell','NaN Matrix','Dates'};
meta.cellType = cellType;
meta.expertise = expertise;
meta.location = location;
meta.depth = depth;



