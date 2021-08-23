% load data
load('Vc.mat','U', 'Vc')
load('opts2.mat','opts')
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask; %maske für cropping

%% alignment zu Allen
U = alignAllenTransIm(U, transParams);
U = rateDisc_removeOutline(U, 10); %ein paar pixels von den seiten weg (grösse sollte 540x640xPCs)
U = arrayCrop(U, allenMask); % trimmt auf die maske 'allenMask' herunter (größe sollte ab hier 540x586xPCs sein)
rateDisc_plotAllenOutline;

%% code to isolate areas
areaIdx = {'VISp' 'SSp-ul' 'MOB' 'SSp-bfd'}; % V1 - HL - M2
targetHS = 'L'; %hemissphere for target area
% get allen area masks
for x = 1:length(areaIdx)
    ind = ismember(dorsalMaps.labelsSplit,areaIdx{x}) & ismember(dorsalMaps.sidesSplit,targetHS);
    areaCoord{x} = poly2mask(dorsalMaps.edgeOutlineSplit{ind}(:,2), dorsalMaps.edgeOutlineSplit{ind}(:,1),size(allenMask,1),size(allenMask,2));
end

%% isolate pixels for target area
x = 4; %show data for V1 (first entry)
areaPixels = arrayShrink(U, ~areaCoord{x}, 'merge'); %this is like U but only for pixels in specific area
areaData = areaPixels * nanmean(Vc,3); %this is all pixels in V1, averaged over all trials
areaPSTH = nanmean(areaData,1); %average over all pixels to get single V1 PSTH
 



