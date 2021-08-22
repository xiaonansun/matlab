function [A,AUC] = twoP_lrCombinAUC_Beta(LR,animal)
% The input LR is a struct array that is the output of
% twoP_lrLoadAllSession.m. The LR array contains multiple fields, of which
% the useful ones include cvAcc, cvAcc_r, mcvAcc_nr_rep, bMaps, and allAUC.
field_name = 'Location';

if ~iscell(animal)
    animal = {animal};
end

T = twoP_getAcquisitionRecord;
iColumn = ismember(T(1,:),{field_name});

A.Setting.location.MM.name = 'MM'; A.Setting.location.MM.iRows = ismember(T(:,iColumn),{A.Setting.location.MM.name});
A.Setting.location.ALM.name = 'ALM'; A.Setting.location.ALM.iRows = ismember(T(:,iColumn),{A.Setting.location.ALM.name});

% Define some settings are constants
A.Setting.baseDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
A.Setting.location.name = 'AllAreas';
A.Setting.depth.sup = [1 224]; % set neuron depth range in microns
A.Setting.depth.int = [225 449]; 
A.Setting.depth.deep = [450 700]; 

A.Setting.animal(1).name = 'CSP27'; A.Setting.animal(1).dateRange = {'3/19/2020','4/1/2020'}; % Define the animal and select a range of dates on which the sessions were performed
A.Setting.animal(2).name = 'CSP30'; A.Setting.animal(2).dateRange = {'3/16/2020','4/1/2020'}; % Define the animal and select a range of dates on which the sessions were performed
A.Setting.animal(3).name = 'Plex50'; A.Setting.animal(3).dateRange = {'3/30/2020','4/2/2020'}; % Define the animal and select a range of dates on which the sessions were performed
A.Setting.animal(4).name = 'Plex51'; A.Setting.animal(4).dateRange = {'3/25/2020','4/2/2020'}; % Define the animal and select a range of dates on which the sessions were performed
A.Setting.animal(5).name = 'Fez51'; A.Setting.animal(5).dateRange = {'4/2/2020','7/16/2020'}; % Define the animal and select a range of dates on which the sessions were performed
A.Setting.animal(6).name = 'Fez57'; A.Setting.animal(6).dateRange = {'6/7/2020','7/16/2020'}; % Define the animal and select a range of dates on which the sessions were performed
A.Setting.animal(7).name = 'Fez59'; A.Setting.animal(7).dateRange = {'6/6/2020','7/3/2020'}; % Define the animal and select a range of dates on which the sessions were performed
iAnimal = find(ismember({A.Setting.animal.name}, animal)); 

A.Date = ['' LR.date];

if isempty(iAnimal)
    disp('This animal does not exist for logistic regression analysis, function is terminating.');
    return;
end

iLR = []; A.currentAnimal = [];
for n = 1:length(iAnimal)
    A.currentAnimal = [A.currentAnimal A.Setting.animal(iAnimal(n)).name];
    idxAnimal = strcmp({LR.animal},A.Setting.animal(iAnimal(n)).name);
    tlower = datetime(A.Setting.animal(iAnimal(n)).dateRange{1},'InputFormat','MM/dd/yyyy');
    tupper = datetime(A.Setting.animal(iAnimal(n)).dateRange{2},'InputFormat','MM/dd/yyyy');
    idxDate = isbetween(A.Date,tlower,tupper);
    A.Setting.animal(iAnimal(n)).iSelect = find(idxAnimal.*idxDate);
    iLR = [iLR A.Setting.animal(iAnimal(n)).iSelect];
end
A.iSelected = iLR;

% A.Setting.animal(iAnimal).dateRange = datetime(A.Setting.animal(iAnimal).dateRange,'InputFormat','MM/dd/yyyy');
% A.Animal = strcmp({LR.animal},A.Setting.animal(iAnimal).name); % Find the index of the rows in the LR structure that corresponds to the specified animal

A.Exist = ~cellfun(@isempty,{LR.cvAcc});
A.Depth = cellfun(@str2double,{LR.depth},'UniformOutput',false);
A.Depth = cell2mat(A.Depth); 
for i = 1:length(A.Depth)
    if isnan(A.Depth(i))
        A.DepthName(i) = {'none'};
    elseif A.Setting.depth.deep(2) >= A.Depth(i) && A.Depth(i) >= A.Setting.depth.deep(1)
        A.DepthName(i) = {'deep'};
    elseif A.Setting.depth.int(2) >= A.Depth(i) && A.Depth(i) >= A.Setting.depth.int(1)
        A.DepthName(i) = {'int'};
    elseif A.Setting.depth.sup(2) >= A.Depth(i) && A.Depth(i) >= A.Setting.depth.sup(1)
        A.DepthName(i) = {'sup'};
    end
end


% A.DateRange = datenum(A.Setting.animal(iAnimal).dateRange(2)) >= datenum(A.Date) & datenum(A.Date) >= datenum(A.Setting.animal(iAnimal).dateRange(1));
iSel = zeros(1,length(LR)); iSel(A.iSelected)=1;
A.iDeep = find(iSel & A.Exist & A.Setting.depth.deep(2) >= A.Depth & A.Depth >= A.Setting.depth.deep(1));
A.iInt = find(iSel & A.Exist & A.Setting.depth.int(2) >= A.Depth & A.Depth >= A.Setting.depth.int(1));
A.iSup = find(iSel & A.Exist & A.Setting.depth.sup(2) >= A.Depth & A.Depth >= A.Setting.depth.sup(1));
A.iMM = find(iSel & A.Exist & A.Setting.location.MM.iRows');
A.iALM = find(iSel & A.Exist & A.Setting.location.ALM.iRows');

% A.iSelected = find(A.Animal & A.DateRange);

A.epoch.iPreStim = [];
A.epoch.iEarlyStim =[];
A.epoch.iLateStim =[];
A.epoch.iDelay =[];
A.epoch.iEarlyRes =[];
A.epoch.iLateRes =[];


%%
AUC.NR.val = []; AUC.R.val = []; 
AUC.NR.valAll = []; AUC.R.valAll = [];
AUC.NR.mu = []; AUC.R.mu = [];
AUC.NR.sigma = []; AUC.R.sigma = [];
AUC.NR.pos = []; AUC.R.pos = [];
AUC.NR.neg = []; AUC.R.neg = [];
AUC.NR.name = {'Non-red'}; AUC.R.name = {'Red'};
AUC.epoch = A.epoch;
AUC.epochNames = fieldnames(A.epoch);

AUC.iSelected = A.iSelected;
AUC.location = A.Setting.location.name;

AUC.MM = AUC; AUC.ALM = AUC.MM;
AUC.MM.iSelected = A.iMM; AUC.ALM.iSelected = A.iALM;
AUC.ALM.location = A.Setting.location.ALM.name; AUC.MM.location = A.Setting.location.MM.name;
% AUC = twoP_computeAUC(AUC,LR);
% AUC.MM = twoP_computeAUC(AUC.MM,LR);
% AUC.ALM = twoP_computeAUC(AUC.ALM,LR);

AUC = twoP_computeAUCFraction(LR,AUC);
AUC.MM = twoP_computeAUCFraction(LR,AUC.MM);
AUC.ALM = twoP_computeAUCFraction(LR,AUC.ALM);

%% Combine AUC and Beta into a single struct
%AUC
A.allAUC.analysisType = 'AUC';

A.allAUC.deep.NR = []; A.allAUC.deep.R = []; A.allAUC.deep.name = 'Deep';
try
for i = 1:length(A.iDeep)
    idx1 = LR(A.iDeep(i)).iEpoch(1,3); idx2 = LR(A.iDeep(i)).iEpoch(2,3);
    A.allAUC.deep.NR = [A.allAUC.deep.NR; mean(LR(A.iDeep(i)).allAUC(LR(A.iDeep(i)).idx_notredcell,idx1:idx2),2)];
    A.allAUC.deep.R = [A.allAUC.deep.R; mean(LR(A.iDeep(i)).allAUC(LR(A.iDeep(i)).idx_redcell,idx1:idx2),2)];
end
end

A.allAUC.int.NR = []; A.allAUC.int.R = []; A.allAUC.int.name = 'Intermediate';
try
for i = 1:length(A.iInt)
    idx1 = LR(A.iInt(i)).iEpoch(1,3); idx2 = LR(A.iInt(i)).iEpoch(2,3);
    A.allAUC.int.NR = [A.allAUC.int.NR; mean(LR(A.iInt(i)).allAUC(LR(A.iInt(i)).idx_notredcell,idx1:idx2),2)];
    A.allAUC.int.R = [A.allAUC.int.R; mean(LR(A.iInt(i)).allAUC(LR(A.iInt(i)).idx_redcell,idx1:idx2),2)];
end
end

A.allAUC.sup.NR = []; A.allAUC.sup.R = []; A.allAUC.sup.name = 'Superficial';
try
for i = 1:length(A.iSup)
    idx1 = LR(A.iSup(i)).iEpoch(1,3); idx2 = LR(A.iSup(i)).iEpoch(2,3);
    A.allAUC.sup.NR = [A.allAUC.sup.NR; mean(LR(A.iSup(i)).allAUC(LR(A.iSup(i)).idx_notredcell,idx1:idx2),2)];
    A.allAUC.sup.R = [A.allAUC.sup.R; mean(LR(A.iSup(i)).allAUC(LR(A.iSup(i)).idx_redcell,idx1:idx2),2)];
end
end
% Beta
A.dBeta.analysisType = 'Beta';

A.dBeta.deep.NR = []; A.dBeta.deep.R = []; A.dBeta.deep.name = 'Deep';
try
for i = 1:length(A.iDeep)
    idx1 = LR(A.iDeep(i)).iEpoch(1,3); idx2 = LR(A.iDeep(i)).iEpoch(2,3);
    A.dBeta.deep.NR = [A.dBeta.deep.NR; mean(LR(A.iDeep(i)).bMaps(LR(A.iDeep(i)).idx_notredcell,idx1:idx2),2)];
    A.dBeta.deep.R = [A.dBeta.deep.R; mean(LR(A.iDeep(i)).bMaps(LR(A.iDeep(i)).idx_redcell,idx1:idx2),2)];
end
end

A.dBeta.int.NR = []; A.dBeta.int.R = []; A.dBeta.int.name = 'Intermediate';
try
for i = 1:length(A.iInt)
    idx1 = LR(A.iInt(i)).iEpoch(1,3); idx2 = LR(A.iInt(i)).iEpoch(2,3);
    A.dBeta.int.NR = [A.dBeta.int.NR; mean(LR(A.iInt(i)).bMaps(LR(A.iInt(i)).idx_notredcell,idx1:idx2),2)];
    A.dBeta.int.R = [A.dBeta.int.R; mean(LR(A.iInt(i)).bMaps(LR(A.iInt(i)).idx_redcell,idx1:idx2),2)];
end
end

A.dBeta.sup.NR = []; A.dBeta.sup.R = []; A.dBeta.sup.name = 'Superficial';
try
for i = 1:length(A.iSup)
    idx1 = LR(A.iSup(i)).iEpoch(1,3); idx2 = LR(A.iSup(i)).iEpoch(2,3);
    A.dBeta.sup.NR = [A.dBeta.sup.NR; mean(LR(A.iSup(i)).bMaps(LR(A.iSup(i)).idx_notredcell,idx1:idx2),2)];
    A.dBeta.sup.R = [A.dBeta.sup.R; mean(LR(A.iSup(i)).bMaps(LR(A.iSup(i)).idx_redcell,idx1:idx2),2)];
end
end