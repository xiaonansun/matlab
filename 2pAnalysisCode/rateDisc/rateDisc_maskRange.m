function [xRange, yRange] = rateDisc_maskRange(mask,minRange)
% get inner range of a binary mask

if ~exist('minRange','var') || isempty(minRange)
    minRange = 0;
end

xMin = find(sum(~mask,1) > 0, 1); %left edge of mask
xMax = find(sum(~mask,1) > 0, 1, 'last'); %right edge of mask
xRange = xMin : xMax; %range to use for X

if length(xRange) < minRange
    missRange = floor((minRange - length(xRange))/2); % get missing range to get to the minimum.
    xRange = xMin-missRange : xMax+missRange;
    if length(xRange) < minRange %for uneven mismatch, reduce by one extra value
        xRange = [xRange xRange(end)+1];
    end
end

yMin = find(sum(~mask,2) > 0, 1); %left edge of mask
yMax = find(sum(~mask,2) > 0, 1, 'last'); %right edge of mask
yRange = yMin : yMax; %range to use for Y

if length(yRange) < minRange
    missRange = floor((minRange - length(yRange))/2); % get missing range to get to the minimum.
    yRange = yMin-missRange : yMax+missRange;
    if length(yRange) < minRange %for uneven mismatch, add one extra value 
        yRange = [yRange yRange(end)+1];
    end
end