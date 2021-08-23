function rateDisc_optoUnilateral(bhv, groupnames)

%% compute detection performance with unilateral stim + delay for each side
fiberColors = {[0 0 1] [1 0 0]};
h1 = figure;
h2 = figure;
Cnt = 0;
for stimLoc = 1 : 2
    for x = 1 : length(bhv)
    
        Cnt = Cnt + 1;
        
        nrMice =  length(unique(bhv{x}.AnimalID)); %number of animals in current group
        cData = NaN(nrMice, 5); %colums 1+2 are for performance with l/r stimulation, 3:5 are contra choice probability
        for iAnimals = 1 : nrMice
            allData = rateDisc_stimANDdelayOptoUnilateral(bhv{x}, bhv{x}.optoPower > 1 & ismember(bhv{x}.AnimalID,iAnimals));
            cData(iAnimals, 1:2) = allData.detect(3) - squeeze(allData.optoDetect(3,stimLoc,:)); %performance change with left or right optogenetic stimulation
            cData(iAnimals, 3:5) = allData.optoContra(stimLoc,:); 
        end
        
        allData = rateDisc_stimANDdelayOptoUnilateral(bhv{x}, bhv{x}.optoPower > 1); %all mice
        
        figure(h1); subplot(2,length(groupnames),Cnt); hold on;
        xx1 = allData.detect(3) - squeeze(allData.optoDetect(3,stimLoc,:)); %optogenetic vs non-optogenetic performance
        xx2 = allData.detect(3) - squeeze(allData.optoDetectUp(3,stimLoc,:)); %same for upper bound
        xx3 = allData.detect(3) - squeeze(allData.optoDetectLow(3,stimLoc,:)); %same for lower bound

        bar(xx1, 'FaceColor', fiberColors{stimLoc});
        errorbar(1:2, xx1, xx2 - xx1, xx1 - xx3, '.k')
        plot(1:2, cData(:,1:2), 'o', 'MarkerFaceColor','w', 'MarkerEdgeColor', 'k'); ylim([-0.1 0.1]);
%         nBarweb(nanmean(cData(:,1:2),1),sem(cData(:,1:2),1),fiberColors{stimLoc}); ylim([0 0.4]);
        
        figure(h2); subplot(2,length(groupnames),Cnt); hold on;
        xx1 = allData.optoContra(stimLoc,1) - 0.5; %contra choices vs chance
        xx2 = allData.optoContraUp(stimLoc,1) - 0.5; %same for upper bound
        xx3 = allData.optoContraLow(stimLoc,1) - 0.5; %same for lower bound
        bar(xx1, 'FaceColor', fiberColors{stimLoc});
        errorbar(1, xx1, xx2 - xx1, xx1 - xx3, '.k');
        
        plot(1, cData(:,3)-0.5, 'o', 'MarkerFaceColor','w', 'MarkerEdgeColor', 'k'); ylim([-0.2 0.1]);
%         %         nBarweb(nanmean(cData(:,3:5),1)-0.5,sem(cData(:,3:5),1),fiberColors{stimLoc}); ylim([-0.1 0.1]);
%         plot(1:3, allData.optoContra(stimLoc,:)-0.5, 'o', 'MarkerFaceColor','w', 'MarkerEdgeColor', 'k');
%         disp(allData.optoContraCnt(stimLoc,1));

    end
end
        