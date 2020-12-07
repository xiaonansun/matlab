%%
clear all;

Animal = 'CSP30';
Session='01-Apr-2020_1';

dataDir='H:\twoP';
lFile={'interpVc'; 'orgdimBeta';'orgfullcorr';'orgregData'};
fileExt='.mat';
for i = 1:length(lFile)
    load([dataDir filesep Animal filesep 'SpatialDisc' filesep Session filesep lFile{i} fileExt]);
end



%%

rho=1;
filterFrames=5;
% regNameIdx=25;
Animal='Plex51';
betaKernel=struct;
% figure(1);
for regNameIdx = 1:length(regLabels)
%     display(['Regressor: ' regLabels{regNameIdx}])
    try
        betaKernel(regNameIdx).value=dimBeta(regIdx(~rejIdx) == regNameIdx,:)';
        betaKernel(regNameIdx).regressorName=regLabels{regNameIdx};
    end
end
betaIdx = find(arrayfun(@(betaKernel) ~isempty(betaKernel.value),betaKernel));


%% Unsorted

figure(1);
suptitle(Animal)
for i = 1:length(betaIdx); 
    regBeta=betaKernel(betaIdx(i)).value;
    [M iM]=max(smoothdata(regBeta,2,'sgolay',filterFrames),[],2);
    % [M iM]=max(betaKernel,[],2);
    iM=sortrows([iM (1:1:length(iM))'],1);
    subplot(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))),i)
    imagesc(regBeta,[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    title(regLabels{betaIdx(i)});
end

%% Sorted
figure(2); 
% suptitle(Animal)
for i = 1:length(betaIdx)
    regBeta=betaKernel(betaIdx(i)).value;
    [M iM]=max(smoothdata(regBeta,2,'sgolay',filterFrames),[],2);
    % [M iM]=max(betaKernel,[],2);
    iM=sortrows([iM (1:1:length(iM))'],1);
    subplot(floor(sqrt(length(betaIdx))),ceil(sqrt(length(betaIdx))),i)
    imagesc(regBeta(iM(:,2),:),[median(regBeta(:))-rho*std(regBeta(:)) median(regBeta(:))+rho*std(regBeta(:))]);
    title(regLabels{betaIdx(i)});
end
