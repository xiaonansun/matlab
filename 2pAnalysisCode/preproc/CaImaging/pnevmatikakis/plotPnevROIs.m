function plotPnevROIs(im, CC, colors)
% plotPnevROIs(backgroundImage, contourOutlines [, colors])
%
% Take the contourOutlines output of ROIContoursPnev, and plots them over
% backgroundImage. If backgroundImage is empty ( [] ), this function will
% plot the outlines into the current figure.
%
% colors is an optional nROIs x 3 array specifying the colors to be used
% for the corresponding ROIs. If not specified or empty, maximally distinct
% colors (from colorcube) will be used.

%% Optional arguments

% Display image if supplied. Otherwise we'll assume this is an existing
% image, and not spawn a new figure
if ~isempty(im)
  figure;
  imagesc(im);
  hold on;
  axis off;
  axis image
%   colormap('bone');
end

% If no colors supplied, use colorcube to get the max number of
% distinguishable colors
if ~exist('colors', 'var') || isempty(colors)
  colors = colorcube(length(CC));
end


%% Plot the ROIs

% Loop through ROIs
for roi = 1:length(CC)
  % Loop through regions of the contours
  % Each contour region starts with a single column of metadata, where the
  % upper element tells you how many points are in the region. (It would
  % normally be the lower element when working with contourc, but we
  % swapped the rows in ROIContoursPnev)
  
  % When a single ROI is composed of two or more separate contour regions,
  % the coordinates of all regions are concatenated in CC, but separated by
  % a column of metadata (FN).
  
  regI = 1; % regI indicate the index of metadata column for each contour region (FN).
  while regI < size(CC{roi}, 2)
    nElem = CC{roi}(1, regI);
    plot(CC{roi}(2, regI + (1:nElem)), CC{roi}(1, regI + (1:nElem)), '-', 'color', colors(roi, :)); % FN removed minus 1 from nElem in order to plot a full circle for the ROI.  
%     plot(CC{roi}(2, regI + (1:nElem-1)), CC{roi}(1, regI + (1:nElem-1)), '-', 'color', colors(roi, :));
    regI = regI + nElem + 1;
  end
end
