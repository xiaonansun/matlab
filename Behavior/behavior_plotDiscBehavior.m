%%
close all;
Animal = 'Plex62';
cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
[Performance,bhv] = DelayedLoc_learningCurves(Animal,cPath);
disp([Animal '. Detection sessions: ' num2str(sum(sum(Performance.Detection>0))) '; Discrimination Sessions: ' num2str(sum(sum(Performance.Discrimination>0)))]);
% [Performance,bhv,allDates] = RateDisc_learningCurves(Animal,cPath);

%% Plot discrimination curves

bhvFileName = [Animal '*Session*.mat'];
bhvFileList = dir(fullfile(cPath,Animal,['**\*' bhvFileName]));
folderList = {bhvFileList.folder}'; filenameList = {bhvFileList.name}';
pathList = cellfun(@(x,y) [x filesep y],folderList,filenameList,'UniformOutput',false);
% d = cellfun(@(x) load(x,'SessionData'),c,'UniformOutput',false);

filenameArray = cell2mat(filenameList);
paradigm = 'SpatialDisc';
cDate = datenum(filenameArray(:,length([Animal '_' paradigm '_'])+1:length([Animal '_' paradigm '_'])+10)); %isolate dates from Filenames
cDate = cDate + (str2num(filenameArray(:,length([Animal '_' paradigm '_'])+19:end-4))*0.01); %add session nr to timestamp
[cDate,ind] = sort(cDate,'descend'); %sort counts to get the order of files to days correct. Newest file should be first in the list.
cDate = floor(cDate);
pathList=pathList(ind,:);
%%
bhv = []; Cnt = 0;
% hF = figure;cl
for iFiles = 1:length(pathList)
   %%
    try
    load(pathList{iFiles}, 'SessionData'); %load current bhv file
    catch
    end
    
    if isfield(SessionData,'Rewarded')
        SessionData.Rewarded = logical(SessionData.Rewarded);
    end
    
    useData = isfield(SessionData,'Rewarded') && length(SessionData.Rewarded) > 100 && sum(SessionData.DistStim) > 0; % use files containing at least 100 trials and different distractor rates
    if useData
        Cnt = Cnt + 1;
        Performance.Date{1,Cnt} = datestr(cDate(iFiles));
        cInd = 1:length(SessionData.Rewarded);
        [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = ...
            rateDisc_audioDiscCurve(SessionData, cInd);
        if isempty(dataUpper) || isempty(dataLower)
            figure
            p = plot(distRatio,pChoseHigh,'.-k');
        else
            [h,p] = boundedline(distRatio,pChoseHigh,[[pChoseHigh - dataLower]; [dataUpper-pChoseHigh]]','.-k','alpha','transparency',0.1);
        end
        
        xlim([0 1]); ylim([0 1]);
        title(['Trials completed: ' num2str(length(SessionData.Rewarded)) '(' num2str(sum(SessionData.Rewarded)) ' rewarded)' ]);
        ylabel('Proportion of right choices');
        xlabel('Right events (%)');
        offsetAxes(gca);
        fig_configAxis(gca);
    end
end

hDisc = figure(3);
cInd = 1:length(SessionData.Rewarded);
[distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh] = ...
    rateDisc_audioDiscCurve(SessionData, cInd);
% rateDisc_audioDiscCurve(bhv, cInd, distBins, discOnly, fixBias, returnCIs)
% plot(distRatio,pChoseHigh,'.-k'); hold on;
% plot(distRatio,dataUpper,'.--k'); 
% plot(distRatio,dataLower,'.--k'); 

[h,p] = boundedline(distRatio,pChoseHigh,[[pChoseHigh-dataLower];[dataUpper-pChoseHigh]]','.-k','alpha','transparency',0.1);

xlim([0 1]); ylim([0 1]);
title([num2str(sum(SessionData.Rewarded)) ' trials rewarded. ' num2str(length(SessionData.Rewarded)) ' trials completed. ']);
ylabel('Proportion of right choices');
xlabel('Right events (%)');
offsetAxes(gca);
fig_configAxis(gca);

%%
cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
animal = 'Plex62';

sDir = dir([cPath filesep animal filesep 'SpatialDisc' filesep 'Session Data' filesep '*.mat']);
for i = 1:length(sDir)
    load([sDir(i).folder filesep sDir(i).name])
    disp([sDir(i).name ': ' num2str(length(SessionData.Rewarded))]);
end