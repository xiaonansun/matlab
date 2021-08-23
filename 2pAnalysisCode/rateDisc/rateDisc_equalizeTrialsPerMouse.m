function ind = rateDisc_equalizeTrialsPerMouse(bhv, ind)
%code to ensure that a similar amount of trials are used from each animal
%in the bhv structure.
% % 
% rng('default');
% tCount = [];
% animals = unique(bhv.AnimalID(ind)); %animals in current selection. equalize trial counts.
% for x = 1:length(animals)
%     tCount(x) = sum(bhv.AnimalID(ind) == animals(x)); %find minimum trialcount
% end
% animals(tCount == 0) = []; %dont consider animals without trials
% for x = 1:length(animals)
%     temp = find(bhv.AnimalID == animals(x) & ind);
%     ind(temp(min(tCount)+1:end)) = false; %reject excess trials
% end