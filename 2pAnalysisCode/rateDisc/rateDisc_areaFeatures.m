function [newX, xArea, xMajor, xMinor, xEccentricity, xOrientation, ...
    xSolidity, xIntensity, xCentroid, xWCentroid] = rateDisc_areaFeatures(X, areaThresh)

%% get image features
xArea = zeros(size(X,3),1,'single');
xMajor = zeros(size(X,3),1,'single');
xMinor = zeros(size(X,3),1,'single');
xEccentricity = zeros(size(X,3),1,'single');
xOrientation = zeros(size(X,3),1,'single');
xSolidity = zeros(size(X,3),1,'single');
xIntensity = zeros(size(X,3),1,'single');
xCentroid = zeros(size(X,3),2,'single');
xWCentroid = zeros(size(X,3),2,'single');
newX = zeros(size(X),'single');

for x = 1 : size(X,3)
    
    rawFrame = X(:,:,x);
    rawFrame(isnan(rawFrame)) = 0;
    
    % smooth and normalize
    cFrame = maxnorm(smoothImg(rawFrame,2,10));
    cFrame = imdilate(cFrame > areaThresh,strel('disk',8));
    cFrame = imerode(cFrame,strel('disk',8));
    cFrame = bwmorph(cFrame,'open'); %break weak connections between different areas
    cFrame = bwareaopen(cFrame,50); %don't use areas that are smaller as 50 pixels
    
    % get area properties
    areaInfo = regionprops(cFrame, rawFrame, ...
        'Centroid', 'WeightedCentroid', 'MeanIntensity', 'Area', 'MajorAxisLength',  ...
        'MinorAxisLength', 'Eccentricity', 'Orientation', 'Solidity', 'Image');
    
    if length(areaInfo) > 1
        areas = cat(1,areaInfo(:).Area); %size of different patches
        cIdx = find((areas - (max(areas) / 2)) > 0); %areas that are at least half as large as the biggest one
        
        if length(cIdx) == 1
            cFrame = bwlabel(cFrame) == cIdx;
            areaInfo = areaInfo(cIdx);
        else
            areaInfo = [];
            cFrame = false(size(cFrame));
        end
    end
    
    rawFrame(~cFrame) = 0;

    % keep variables
    if ~isempty(areaInfo)
        xArea(x) = areaInfo.Area;
        xMajor(x) = areaInfo.MajorAxisLength;
        xMinor(x) = areaInfo.MinorAxisLength;
        xEccentricity(x) = areaInfo.Eccentricity;
        xOrientation(x) = areaInfo.Orientation;
        xSolidity(x) = areaInfo.Solidity;
        xIntensity(x) = areaInfo.MeanIntensity;
        xCentroid(x,:) = areaInfo.Centroid;
        xWCentroid(x,:) = areaInfo.WeightedCentroid;
        newX(:,:,x) = rawFrame; %keep thresholded area
    end
end
end