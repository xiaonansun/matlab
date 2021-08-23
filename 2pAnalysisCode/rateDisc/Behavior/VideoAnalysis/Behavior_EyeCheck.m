function [params, didFuse] = Behavior_EyeCheck(cFrame,highThresh, lowThresh, showPlot, forceFuse, fusePower)

if ~exist('lowThresh', 'var') || isempty(lowThresh)
    lowThresh = 0;
end

if ~exist('showPlot', 'var') || isempty(showPlot)
    showPlot = false;
end

if ~exist('forceFuse', 'var') || isempty(forceFuse)
    forceFuse = false;
end

if ~exist('fusePower', 'var') || isempty(fusePower)
    fusePower = 2;
end

if showPlot
    imshow(cFrame);
end

didFuse = false; %flag for scattered image. This is true when code tried to fuse small pieces further togegher.

%% improve image contrast, remove noise and close gaps
% cFrame = imadjust(cFrame,[lowThresh graythresh(cFrame)/3],[]); %improve contrast
cFrame = imadjust(cFrame,[lowThresh highThresh],[]); %improve contrast
level = graythresh(cFrame);
cFrame = cFrame < level;
cFrame = bwareaopen(cFrame,10); %don't use areas that are smaller as 10 pixels
areaInfo = regionprops(cFrame, 'Area');

if forceFuse
    if fusePower > 0 %this will fuse patches
        cFrame = imdilate(cFrame,strel('disk',fusePower));
        cFrame = imerode(cFrame,strel('disk',fusePower));
    elseif fusePower < 0 %this will erode patches more strongly
        cFrame = imerode(cFrame,strel('disk',-fusePower));
        cFrame = imdilate(cFrame,strel('disk',-fusePower));
    end
    didFuse = true;
else
    if (length(areaInfo) > 3 && ~any([areaInfo.Area] > 150)) %image is very scattered. Allow some areas to fuse together.
        cFrame = imdilate(cFrame,strel('disk',2));
        cFrame = imerode(cFrame,strel('disk',2));
        didFuse = true;
    elseif any([areaInfo.Area] > 700) %large area is present do some erosion to make sure this is not a false merge.
        cFrame = imerode(cFrame,strel('disk',4));
        cFrame = imdilate(cFrame,strel('disk',4));
        didFuse = true;
    end
end

cFrame = bwmorph(cFrame,'open'); %break weak connections between different areas
cFrame = imclearborder(cFrame); %don't use areas that touch image border
areaInfo = regionprops(cFrame, 'Centroid');

if length(areaInfo) > 1 %if more than one area is present, find combination that is closest to circular 
    cntr = .5 * [size(cFrame,2) size(cFrame,1)]; % X-Y coordinates and NOT Row/Col
    d = sqrt( sum( bsxfun(@minus,vertcat( areaInfo.Centroid ), cntr ).^2, 2 ) );
    [~, idx] = min(d); %index for area that is closest to the center
    checkFrame = bwlabel(cFrame);
    tIdx = 1:length(areaInfo); %all areas
    tIdx = [idx tIdx((tIdx ~= idx))]; %center area should come first

    for iCombs = 1:length(areaInfo) + 1
        if iCombs == length(areaInfo) + 1
            temp = checkFrame > 0;
        else
            temp = checkFrame == idx | checkFrame == tIdx(iCombs); %combine center with other areas to check if it gives better result
        end
        
        for iTest = 1:50 %make sure areas get merged properly
            areaTest = regionprops(imclose(temp,strel('disk',20+iTest)), 'Eccentricity', 'Area','Solidity');
            if length(areaTest) == 1
                break
            end
        end
        allEcc(iCombs) = areaTest.Eccentricity; %eccentricity is the measure for roundness here. lower values are better.
        allArea(iCombs) = areaTest.Area; %size of current fit.
        allSolid(iCombs) = areaTest.Solidity; %Solidity of current fit.
    end
       
    [minEcc, idx] = min(allEcc); %index for combination that ends up most circular
    
    if any(allSolid(allEcc > minEcc) > 0.85) && allArea(idx) < 100 %if other decent options are available and current selection is very small, use next better option instead
        allEcc(idx) = 1;
        [~, idx] = min(allEcc); %next best index
    end    
    
    if idx ~= 1
        if idx < length(allEcc)
            idx = [1 idx];
        else
            idx = 1:length(tIdx);
        end
    end
    cFrame = ismember(checkFrame,tIdx(idx)); %return best combination and continue
end

areaInfo = regionprops(cFrame, 'Area', 'Centroid', 'MajorAxisLength','MinorAxisLength','Orientation','Solidity');
if length(areaInfo) > 1
    for iTest = 1:50
        areaInfo = regionprops(imclose(cFrame,strel('disk',20+iTest)), 'Area', 'Centroid', 'MajorAxisLength','MinorAxisLength','Orientation','Solidity');
        
        if length(areaInfo) == 1
            break
        end
    end
end

if length(areaInfo) == 1
    params.axes = [areaInfo.MajorAxisLength areaInfo.MinorAxisLength];
    params.orient = areaInfo.Orientation;
    params.solidity = areaInfo.Solidity;
    params.center = areaInfo.Centroid;
else
    params.orient = 0;
    params.axes = zeros(1,2);
    params.center = zeros(1,2);
    params.solidity = 0;
end
 
%% find pupil ellipse
if showPlot && sum(params.center(:)) ~= 0
    phi = linspace(0,2*pi,50);
    cosphi = cos(phi);
    sinphi = sin(phi);
    
    theta = pi*params.orient/180;
    R = [ cos(theta)   sin(theta)
        -sin(theta)   cos(theta)];
    
    xy = [(params.axes(1)/2)*cosphi; (params.axes(2)/2)*sinphi];
    xy = R*xy;
    
    x = xy(1,:) + params.center(1);
    y = xy(2,:) + params.center(2);

    hold on
    plot(x,y,'r','LineWidth',2);
    plot(params.center(1),params.center(2),'xr','LineWidth',2);
    hold off
end
end

% %% nested code
% function params = getPupil(frameIn)
% areaInfo = regionprops( frameIn, 'Area', 'Centroid', 'MajorAxisLength','MinorAxisLength','Orientation','Solidity');
% 
% %find region closest to center
% if ~isempty(areaInfo)
%     cntr = .5 * [size(frameIn,2) size(frameIn,1)]; % X-Y coordinates and NOT Row/Col
%     d = sqrt( sum( bsxfun(@minus,vertcat( areaInfo.Centroid ), cntr ).^2, 2 ) );
%     [~, idx] = min(d);
%     
%     params.axes = [areaInfo(idx).MajorAxisLength areaInfo(idx).MinorAxisLength];
%     params.orient = areaInfo(idx).Orientation;
%     params.solidity = areaInfo(idx).Solidity;
%     params.center = areaInfo(idx).Centroid;
% else
%     params.orient = 0;
%     params.axes = zeros(1,2);
%     params.center = zeros(1,2);
%     params.solidity = 0;
% end
% end
    
%% old code
%disable warnings from imfindcircles
% warning('off','images:imfindcircles:warnForSmallRadius');
% warning('off','images:imfindcircles:warnForLargeRadiusRange');

% cFrame = imadjust(cFrame,[0 imThresh],[]); %improve contrast
% % cFrame = conv2(double(cFrame), Kernel, 'valid'); %convolve original matrix with filter
% level = graythresh(cFrame);
% cFrame = cFrame < level;
% cFrame = bwareaopen(cFrame,10);
% cFrame = imclose(cFrame,strel('disk',2));
% cFrame = imerode(cFrame,strel('disk',2));
% cFrame = imdilate(cFrame,strel('disk',2));
%
% [params.center, params.radius] = imfindcircles(cFrame,[3 20]);
% %find region closest to center
% if ~isempty(params.center)
%     cntr = .5 * [size(cFrame,2) size(cFrame,1)]; % X-Y coordinates and NOT Row/Col
%     d = sqrt( sum( bsxfun(@minus,vertcat( params.center ), cntr ).^2, 2 ) );
%     [~, idx] = min(d);
%
%     params.center = params.center(idx,:);
%     params.radius = params.radius(idx);
% else
%     params.center = zeros(1,2);
%     params.radius = 0;
% end
%
% % [params.center, params.radius] = imfindcircles(cFrame,[max([round(params.radius)-5 3]) round(params.radius)+5]);
% % %find region closest to center
% % if ~isempty(params.center)
% %     cntr = .5 * [size(cFrame,2) size(cFrame,1)]; % X-Y coordinates and NOT Row/Col
% %     d = sqrt( sum( bsxfun(@minus,vertcat( params.center ), cntr ).^2, 2 ) );
% %     [~, idx] = min(d);
% %
% %     params.center = params.center(idx,:);
% %     params.radius = params.radius(idx);
% % else
% %     params.center = zeros(1,2);
% %     params.radius = 0;
% % end
%
% if showPlot
% viscircles(params.center, params.radius,'LineStyle','--','LineWidth',1);
% end