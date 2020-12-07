function rateDisc_plotAreaSVD(animal)
%% some variables
% nDims = [100, 500, 1000]; %number of whole-frame components to test
if ispc
    cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
dPath = [cPath animal filesep 'SpatialDisc' filesep]; % path for raw data
noiseThresh = 0.5; %threshold for mse to reject bad recording

%% load explained ariance results
errFiles = dir([bPath 'mse_*']);
for iFiles = 1 : length(errFiles)
    temp = textscan(errFiles(iFiles).name, '%s%d%s', 'Delimiter','_');
    nDims(iFiles) = temp{2};
    
    load([errFiles(iFiles).folder filesep errFiles(iFiles).name])
    allMse(iFiles,:) = mse;
end

%% plot overview
cMap = winter(length(allMse)+1);
useIdx = allMse(nDims == 200, :) < noiseThresh;
fprintf('Rejected %d/%d recordings for mse>%g\n',sum(~useIdx),length(useIdx),noiseThresh);
allMse(:,~useIdx) = [];

[nDims,b] = sort(nDims);
allMse = allMse(b,:);
figure; 
subplot(1,2,1);hold on
plot(allMse(nDims == 200,:)); axis square; ylabel('mse'); xlabel('recordings'); title(['Prediction error; 200PCs - ' animal])
xlim([0 sum(useIdx)+1]);
subplot(1,2,2);hold on
for iFiles = 1 : length(allMse)
    plot(nDims, (1-allMse(:,iFiles)),'o-', 'linewidth',2, 'color',cMap(iFiles+1,:));
    ylabel('explained variance'); xlabel('#PCs'); title(['Dimensionality - ' animal]);
    xlim([0 nDims(end)+nDims(1)]);
end
plot(nDims, (mean(1-allMse,2)),'ko-', 'linewidth',4,'MarkerSize',10); hline(0.9);
axis square

end


