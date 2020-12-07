%% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign,cPath] = twoP_delayDecRecordings;
taskLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

%% get reconstructed Vs, modality index, datapath and allOpts
[recV, semV, allBeta, recLabels, dataPath, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes, cellCnt] = twoP_motorReconstruct('All'); %get reconstructed V, used full model

recV = cat(4,recV{:});
stimTimes = cat(2,stimTimes{:});
stimTimes(stimTimes > 1.5) = [];
figure; histogram(stimTimes,50); axis square; %show stim time histogram
nRecLabels = [{'full'} recLabels{:} {'motor' 'motorVideo' 'sensory' 'cognitive'}];
xRange = [(-baseLength : 1 : 0)./31, (1 : frames-baseLength-1) ./31]; %time vector for x-axis. Stimulus onset is at 0s.

%% compute modulation index
fullRec = sum(recV,4);
motorRec = sum(recV(:,:,:,ismember(recLabels,motorLabels)),4);
taskRec = sum(recV(:,:,:,ismember(recLabels,taskLabels)),4);

iMod = 6; %all trials
allVar = sum(abs(fullRec(:,:,iMod)),2);
motorVar = sum(abs(motorRec(:,:,iMod)),2);
taskVar = sum(abs(taskRec(:,:,iMod)),2);

modIndex = ((taskVar-motorVar) ./ (taskVar + motorVar) + 1) / 2;

% make figure
[a, b] = (sort(modIndex,'ascend')); %sort indices
ind = round(a*200);
modMat = false(size(modIndex,1),200);

for x = 1 : length(ind)
    modMat(x, 1:ind(x)) = true;
end

figure
imagesc((modMat));
caxis([-0.5 1.5]);
vline([mean(ind) 100],'w');
colormap(colormap_blueblackred);
axis square

figure;
histogram(modIndex,0:0.05:1)
vline(median(modIndex),'r'); axis square
title(['Median modIdx: ' num2str(median(modIndex))])
    
    
%% trace examples
[~,c] = sort(allVar,'descend');
figure
for x = 1:100
plot(fullRec(c(x),:,6))
hold on; axis square
plot(motorRec(c(x),:,6))
plot(taskRec(c(x),:,6))
vline([61 81 81+round((1.7*31)) 81+round((2.7*31))],'k')
hold off
pause
end