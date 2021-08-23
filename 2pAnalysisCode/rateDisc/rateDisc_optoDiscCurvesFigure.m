function rateDisc_optoDiscCurvesFigure(bhv,groupnames,distBins)
%code to to fit sigmoids to discrimination data and check for the impact of
%optogenetic inhbition.

h = figure;
fiberColors = {[0 0 1] [1 0 0]};
ctrlColor = zeros(1,3);
fiberLocation = {'Frontal' 'Parietal'};

%% discrimination curves
Cnt = 0;
for cLoc = 1 : 2
    for x = 1 : length(bhv)
        Cnt = Cnt + 1;
        clear cData
        oInd = (bhv{x}.optoDur == 1.3 | bhv{x}.optoDur == 1.5 | bhv{x}.optoDur == 0.6) & bhv{x}.optoSide == 3 & bhv{x}.optoPower > 1; %bilateral stim and delay stimuation at max. power.
        pInd = ~bhv{x}.DidNotChoose & ~bhv{x}.DidNotLever & logical(bhv{x}.Assisted & bhv{x}.DistStim > 0); %only use active discrimination trials
        
        clear cInd
        cInd(1,:) = pInd & oInd & bhv{x}.stimLocation == cLoc & bhv{x}.optoType == 1; %optogenetic; frontal stim only
        cInd(2,:) = pInd & bhv{x}.optoDur == 0 & ismember(bhv{x}.SessionNr, unique(bhv{x}.SessionNr(cInd(1,:)))); %non-optogenetic trials from same sessions

        % get absolute performance
        [tOptError{cLoc}(x,1), tOptError{cLoc}(x,2), tOptDetect{cLoc}(x)] = Behavior_wilsonError(sum(bhv{x}.Rewarded & cInd(1,:)), sum(cInd(1,:))); %opto performance
        [tCtrlError{cLoc}(x,1), tCtrlError{cLoc}(x,2), tCtrlDetect{cLoc}(x)] = Behavior_wilsonError(sum(bhv{x}.Rewarded & cInd(2,:)), sum(cInd(2,:))); %control performance
        
        % plot fit and real data
        subplot(2,length(bhv),Cnt); hold on
        [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh, stats] = rateDisc_audioDiscCurve(bhv{x}, cInd(2,:), distBins, false, true, true);
        ctrlParams{cLoc}(x,:) = struct2array(params); %fitted paramters
        ctrlParamsBstr{cLoc}(x,:,:) = struct2array(stats.paramBootstraps); %bootstrap results
        ctrlParamsSE{cLoc}(x,:) = struct2array(stats.paramCIs) ./ 1.96 ; %confidence intervals
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', ctrlColor);
        h(x) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', ctrlColor, 'MarkerFaceColor','w','linewidth',2);
        
        [distRatio, rightChoice, nTrials, params, cFit, dataUpper, dataLower, pChoseHigh, stats] = rateDisc_audioDiscCurve(bhv{x}, cInd(1,:), distBins, false, true, true);
        optParams{cLoc}(x,:) = struct2array(params); %fitted paramters
        optParamsBstr{cLoc}(x,:,:) = struct2array(stats.paramBootstraps); %bootstrap results
        optParamsSE{cLoc}(x,:) = struct2array(stats.paramCIs) ./ 1.96 ; %confidence intervals
        plot(cFit(1,:),cFit(2,:),'linewidth', 4, 'color', fiberColors{cLoc});
        h(x) = errorbar(distRatio, pChoseHigh, pChoseHigh-dataLower, dataUpper-pChoseHigh, 'o', 'color', fiberColors{cLoc}, 'MarkerFaceColor','w','linewidth',2);
        
        xlim([min(distRatio)-mean(diff(distRatio)),max(distRatio)+mean(diff(distRatio))]);
        ylim([0 1]); ylabel('Proportion chose right');axis square
        xlabel('Distractor ratio'); title(groupnames{x});
        
        % add lines for individual mice and keep params
        nCnt = 0;
        for iAnimals = unique(bhv{x}.AnimalID)
            aInd = ismember(bhv{x}.AnimalID, iAnimals);
            
            tInd = false(1,length(bhv{x}.Rewarded));
            temp = find(cInd(1,:));
            tInd(temp(randperm(length(temp),round(length(temp)*0.75)))) = true;
            [distRatio, ~, nTrials, oParams, cFit,~,~,pChoseHigh] = rateDisc_audioDiscCurve(bhv{x}, tInd & aInd, distBins, false, true);
            

            [distRatio, ~, nTrials, oParams, cFit,~,~,pChoseHigh] = rateDisc_audioDiscCurve(bhv{x}, cInd(1,:) & aInd, distBins, false, true);
            if ~isempty(cFit)
                nCnt = nCnt + 1;
                paramNames = fieldnames(oParams);
                plot(cFit(1,:),cFit(2,:),'linewidth', 2, 'color', [fiberColors{cLoc} 0.25]);
                [~, ~, ~, cParams, cFit,~,~,pChoseHigh] = rateDisc_audioDiscCurve(bhv{x}, cInd(2,:) & aInd, distBins);
                plot(cFit(1,:),cFit(2,:),'linewidth', 2, 'color', [ctrlColor 0.5]);
                
                [discOptError{cLoc,x}(nCnt,1), discOptError{cLoc,x}(nCnt,2), discOptPerf{cLoc,x}(nCnt)] = Behavior_wilsonError(sum(bhv{x}.Rewarded & cInd(1,:) & aInd), sum(cInd(1,:) & aInd)); %opto performance
                [discCtrlError{cLoc,x}(nCnt,1), discCtrlError{cLoc,x}(nCnt,2), discCtrlPerf{cLoc,x}(nCnt)] = Behavior_wilsonError(sum(bhv{x}.Rewarded & cInd(2,:) & aInd), sum(cInd(2,:) & aInd)); %control performance
                
            else
                fprintf('No %s discrimination from %s\n',fiberLocation{cLoc},bhv{x}.Animals{iAnimals});
            end
        end
    end
end

%% compare model parameters for lapses and sensitivity
figure;
Cnt = 0;
for cLoc = 1 : 2
    for x = 1 : 3
        Cnt = Cnt +1;
        ax = subplot(2,3,Cnt); hold on;
        
        cData = [ctrlParams{cLoc}(:,x+1),optParams{cLoc}(:,x+1)]; %combine ctrl and test parameter results
        cError = [ctrlParamsSE{cLoc}(:,x+1),optParamsSE{cLoc}(:,x+1)]; %combine ctrl and test parameter results
        nBarweb(cData, cError); %show error bar
        
        ax.XTick = 1 : length(groupnames);
        ax.XTickLabel = groupnames; axis square; xlim([0.5 4.5])
        ax.TickLength = [0 0]; title(paramNames{x+1});
        
    end
end

%% compute total change in performance by collapsing over all distractor rates
figure;
for cLoc = 1 : 2
    ax = subplot(1,2,cLoc); hold on;
    
    Cnt = 1;
    for x = 1 : length(bhv)
        cData = (discCtrlPerf{cLoc,x}' - discOptPerf{cLoc,x}')*100; %combine ctrl and test parameter results
        errorbar(Cnt,nanmean(cData),sem(cData),'k','linewidth',2);
        bar(Cnt,nanmean(cData));
        plot(Cnt,cData,'ko','MarkerFaceColor','w')
        Cnt = Cnt + 1;
    end
    
    ax.XTick = 1 : length(groupnames);
    ax.XTickLabel = groupnames; axis square; ylim([-2 20])
    ax.TickLength = [0 0]; title(fiberLocation{cLoc});
    ylabel('Performance reduction(%)');
    niceFigure;
end
