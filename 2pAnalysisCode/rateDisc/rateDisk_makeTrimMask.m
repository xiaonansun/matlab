Cnt = 0;
nAreaMap = areaMap;

% merge some areas together
nAreaMap(ismember(nAreaMap(:),[51 53 55])) = 51; %merge left retrosplinal cortex
nAreaMap(ismember(nAreaMap(:),[52 54 56])) = 52; %merge right retrosplinal cortex
nAreaMap(ismember(nAreaMap(:),[47 49])) = 47; %merge left AAC
nAreaMap(ismember(nAreaMap(:),[48 50])) = 48; %merge right AAC
nAreaMap(ismember(nAreaMap(:),[23 25 27 29])) = 23; %merge left auditory cortices
nAreaMap(ismember(nAreaMap(:),[24 26 28 30])) = 24; %merge right auditory cortices
nAreaMap(ismember(nAreaMap(:),[21 23 25 27 29])) = 23; %merge left auditory cortices
nAreaMap(ismember(nAreaMap(:),[22 24 26 28 30])) = 24; %merge right auditory cortices
nAreaMap(ismember(nAreaMap(:),[1 5])) = 5; %merge frontal cortices
nAreaMap(ismember(nAreaMap(:),[2 6])) = 6; %merge frontal cortices

trimMap = zeros(size(nAreaMap));
for iAreas = unique(nAreaMap(~isnan(nAreaMap)))'
    Cnt=  Cnt +1;

    cFrame = nAreaMap == iAreas;
    cFrame = imclose(cFrame,strel('disk',4));
    cFrame = bwareaopen(cFrame,100); %don't use areas that are smaller as 100 pixels
    trimMap(cFrame) = iAreas;
    
end


%% make regionmap
load C:\Users\smusall\Documents\repoland\playgrounds\Simon\trimMap.mat
trimMap = trimMap(1:540,:);
regions = {[63 64], [3:6 47 48], 7:20, [23 24 43:46 61 62], [31:42 51 52 57:60]};

regionMap = zeros(size(trimMap));
for x = 1 : length(regions)
    cFrame = ismember(trimMap, regions{x});
    cFrame = imclose(cFrame,strel('disk',4));
    contour(cFrame, 'linewidth',4);
    regionMap(cFrame) = x;
end
    




