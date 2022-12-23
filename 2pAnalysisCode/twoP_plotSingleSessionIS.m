function twoP_plotSingleSessionIS(animal,session)
%%

% Constants
% animal = 'CSP30'; session = '200330';
colAnimal = 1;
colLocation = 3;
colDepth = 4;
colSession = 6;

exps = twoP_getAcquisitionRecord; % load acquisition record
idxRow = contains(exps(:,colAnimal),animal) & contains(exps(:,colSession),session);
location = exps{idxRow,colLocation};
depth = exps{idxRow,colDepth};

S = twoP_settings;
imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;

sPerFrame = S.msPerFrame/1000;
sides = {'Left';'Right'};

% Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'Vc.mrat'),'Vc'); Vc = Vc.Vc;
Vc = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'data.mat'),'data'); Vc = Vc.data.sdu;
cBhv = load(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat'),'cBhv'); cBhv = cBhv.cBhv;
idxCell = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'iscell.npy'));
idxRed = readNPY(fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'redcell.npy'));
idxRed = idxRed(find(idxCell(:,1)));

sBhv = twoP_bhvSubSelection(cBhv); %% Load trial-matched Bpod behavior data

allNan = sum(sum(isnan(Vc),3)); 
filterNan = ones(1,length(allNan)); 
[valPks,idxPks] = findpeaks(allNan);
mu = 0.75;
for i = 1:length(idxPks)
    [valMin(i),idxMin(i)]=min(abs(allNan(S.segFrames(i):idxPks(i))-mu*valPks(i)));
    filterNan(S.segFrames(i)+idxMin(i):idxPks(i))=nan;
end
fVc = Vc.*filterNan; %%

t = 0:sPerFrame:sPerFrame*(size(fVc,2)-1);
transparency = 0.1;
clrOrange = [0.9290 0.6940 0.1250]; clrBlue = [0 0 1]; clrCyan = [0 1 1]; clrMagenta = [1 0 1];
lineColor = [clrOrange;clrMagenta];
idxLeftTrials = contains(sBhv.sub.AllNames,'Left');
idxLeftTrials = find(idxLeftTrials);
nTrialTypes = length(idxLeftTrials);

figAllTrialTypes= figure(1);
set(figAllTrialTypes,'Position',[50 50 500 800]);

tileAllTrialTypes = tiledlayout(nTrialTypes,2, ...
    'TileIndexing','rowmajor',...
    'TileSpacing','tight',...
    'Padding','tight');
idxTrials = zeros(nTrialTypes,2);
for i = 1:nTrialTypes % stim vs response vs reward vs error trials: rows of the tiled plot
    for k = 1:2 % red versus non-red cells: columns of the tiles
        nexttile
        lgdLines = [];
        for j = 1:length(sides) %left versus right: the two traces of a single tile
            if j == 1
                idxTrials(i,j) = idxLeftTrials(i);
            elseif j == 2
                idxTrials(i,j) = idxLeftTrials(i)+1;
            end
            if k == 1
                fVcMean = mean(fVc(find(idxRed),:,sBhv.sub.AllIdx(idxTrials(i,j),:)'),3,'omitnan');
                [lineOut{j},~] = stdshade(fVcMean,transparency,lineColor(j,:),t,[]);
            elseif k == 2
                fVcMean = mean(fVc(find(~idxRed),:,sBhv.sub.AllIdx(idxTrials(i,j),:)'),3,'omitnan');
                [lineOut{j},~] = stdshade(fVcMean,transparency,lineColor(j,:),t,[]);
            end
            iTileTitle = strfind(sBhv.sub.AllNames{idxTrials(i,j)},sides{j});
            tileTitle = sBhv.sub.AllNames{idxTrials(i,j)}(1:iTileTitle-1);
            title(tileTitle)
            ax = gca;
            ax = fig_configAxis(ax);
            lgdLines = [lgdLines lineOut{j}];
        end
        lgd = legend(lgdLines,{sides{1};sides{2}},...
            'Location','Northwest',...
            'Box','off');
        if i ~= nTrialTypes; ax.XTickLabel = []; ax.XAxis.Color = 'none'; end
        
    end
    xlim([min(t) max(t)]);
end
title(tileAllTrialTypes,[animal ' ' session ' ' location ' ' depth '\mum']);
xlabel(tileAllTrialTypes,'Time (s)');
ylabel(tileAllTrialTypes,{'Inferred spikes'});
figureSavePath = fullfile(S.dir.imagingRootDir,'PETH',[animal '_' session '_IS.pdf']);
exportgraphics(figAllTrialTypes,figureSavePath);
disp(['Figure saved as: ' figureSavePath]);
