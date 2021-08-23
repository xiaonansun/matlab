function [clusters, model, objs, colors] = rateDisc_clustersFromGMM(pts, model, useLABColor)
% [clusters, model, objs, colors] = clustersFromGMM(pts)
% [clusters, model, objs, colors] = clustersFromGMM(pts, model [, useLABColor])
% [clusters, model, objs, colors] = clustersFromGMM(pts, kToTry [, useLABColor])
% 
% Take the points with coordinates pts (nPts x nDims), run DBSCAN to get
% rid of unclusterable points, then fit a Gaussian Mixture Model. Sweep the
% number of clusters and use BIC to determine the "correct" number of
% clusters.
% 
% DBSCAN is run with epsilon = range * 0.05. This won't necessarily perform
% well for identifying clusters, but it doesn't matter because all we need
% is for it to find the dense regions so we can run the GMM.
% 
% Outputs are the cluster identities (nPts x 1), best GMM model, the cell
% array of all GMM models, and "colors". "model" may be fed back to get the
% plots and skip the fitting. "colors" is found based on the cluster
% centroids superimposed on the HSV color wheel, and can be fed into
% plotCaMDS in lieu of clusters.
% 
% If kToTry is a vector of numbers, these numbers of clusters will be
% tried. If empty, tests 2:12 clusters. If model is supplied, skips fitting
% and just makes the plots. Points more than 4 SDs away from any cluster
% are also considered background.
% 
% Makes 3 plots: (1) a DBSCAN plot, where non-background points are black
% and background points are gray; (2) contours of the Gaussians over the
% points; (3) the points colored by group membership.


%% Parameters

maxClusters = 50;   % Max number of clusters to try
reps = 10;          % Number of times to try fitting each GMM model
maxIter = 1000;     % Max number of steps to get GMM to converge
whichIC = 'BIC';    % Information criterion to use. BIC works much better than AIC (which doesn't penalize enough)
DBSCANEps = 0.07;   % Scaling factor relative to range for finding dense regions
DBSCANMinPts = 20;  % Minimum number of points per DBSCAN cluster. Must be bigger than the little clumps.
maxMahal = 4;       % If a point is farther than this many SDs from any GMM cluster, it's discarded
contourRes = 100;   % For drawing Gaussians figure

% Can set this to 1 to color clusters based on their centroids' positions
% on the color wheel (to make clear similar vs. distant clusters)
cmapByCentroid = 1;

% These extra parameters are only needed when cmapByCentroid == 1
colorRangeQuantile = 0.01;    % For centering the low-D plot on the color wheel
colorScaleQuantile = 0.98;    % For discarding outliers when imposing the low-D on color wheel

% LAB color related parameters
if ~exist('useLABColor', 'var')
  useLABColor = 0;
end
showLabSpace = 1;
% colorwheelRes is also used for estimating the LAB space borders
colorwheelRes = 400;
LABTol = 1e-3;
LLevel = 60;
chromaRange = 53;



warning('Possible bug: dbIsNoise isn''t passed back; does these mean that computing vs. replotting a model produce slightly different results?');



%% Optional argument

if ~exist('model', 'var')
  model = [];
end


if isempty(model) || isnumeric(model)
  %% Use DBSCAN if requested
  
  if DBSCANEps > 0
    ranges = range(pts);
    
    [~, dbIsNoise] = DBSCAN(pts, DBSCANEps * max(ranges), DBSCANMinPts);
    
    % Plot used points in black, unused points in gray
    figure;
    hold on;
    plot(pts(dbIsNoise, 1), pts(dbIsNoise, 2),  'o', 'MarkerFaceColor', 0.7 * [1 1 1], ...
      'MarkerSize', 6, 'MarkerEdgeColor', 'none');
    plot(pts(~dbIsNoise, 1), pts(~dbIsNoise, 2),  'o', 'MarkerFaceColor', 'k', ...
      'MarkerSize', 6, 'MarkerEdgeColor', 'none');
    set(gca, 'TickDir', 'out', 'Box', 'off');
    axis tight;
    blankAxes;
    
  else
    dbIsNoise = false(1, size(pts, 1));
  end
  
  
  %% Fit the Gaussian Mixture models
  
  options = statset('MaxIter', maxIter);
  
  objs = cell(1, maxClusters);
  
  fitPts = pts(~dbIsNoise, :);
  
  if isempty(model)
    kToTry = 1:maxClusters;
  else
    kToTry = model;
  end
  
  for k = kToTry
    objs{k} = gmdistribution.fit(fitPts, k, 'Replicates', reps, 'Options', options);
  end
  
  % The below code, if uncommented, seeds the fitting algorithm with a
  % guess that has a very broad "background" component to try to suck up
  % the un-clusterable points. In practice this doesn't seem to work.
%   guessSigma = 0.2 * range(pts(:, 1));
%   for k = 1:maxClusters
%     NlogL = Inf;
%     
%     for rep = 1:reps
%       % Construct starting guess
%       thesePts = randperm(size(pts, 1), k);
%       guess.mu = pts(thesePts, :);
%       guess.Sigma = guessSigma * repmat(eye(2), [1 1 k]);
%       guess.Sigma = bsxfun(@times, guess.Sigma, 0.01 * range(pts(:, 1)) + rand([1 1 k]));
%       guess.PComponents = ones(1, k);
%       
%       % Make first component the background uniform-ish component
%       guess.mu(1, :) = [0 0];
%       guess.Sigma(:, :, 1) = 1000 * eye(2);
%       
%       obj = gmdistribution.fit(pts, k, 'Start', guess, 'Options', options);
%       
%       if obj.NlogL < NlogL
%         objs{k} = obj;
%         NlogL = obj.NlogL;
%       end
%     end
%   end
  
  
  %% Find best number of components using BIC
  
  IC = cellfun(@(m) m.(whichIC), objs);
  
  [~, nClusters] = min(IC);
  
  figure;
  plot(IC);
  xlabel('# components');
  ylabel(whichIC);
  set(gca, 'TickDir', 'out', 'Box', 'off');
  
  model = objs{nClusters};

else
  objs = [];
end


nClusters = model.NComponents;

%% Repeat clustering using best model

[clusters, ~, ~, ~, mahal] = cluster(model, pts);


%% Plot points over contours

mins = min(pts);
maxes = max(pts);

blankFigure;
hold on;
plot(pts(:, 1), pts(:, 2),  'o', 'MarkerFaceColor', 'k', ...
  'MarkerSize', 6, 'MarkerEdgeColor', 'none');

[X, Y] = meshgrid(linspace(mins(1), maxes(1), contourRes), linspace(mins(2), maxes(2), contourRes));
Z = reshape(pdf(model, [X(:), Y(:)]), size(X));
contour(X, Y, Z, 30);

axis tight;
axis square;


%% Clear classification of uncertain points

clusters(min(mahal, [], 2) > maxMahal) = 0;


%% Prep for using LAB color

if useLABColor
  % Create an oversize LAB space image
  overChromaRange = 100;
  wheelImLab = createLABIm(LLevel, overChromaRange, colorwheelRes);
  
  wheelIm = mapLABToRGBAndBack(wheelImLab, LABTol);
  
  % Find COM
  [i1, i2] = ind2sub(size(wheelIm), find(~isnan(wheelIm(:, :, 1))));
  imCSCOM(1) = mean(i1);
  imCSCOM(2) = mean(i2);
  
  % Recenter, rescale COM
  CSCOM = imCSCOM - colorwheelRes / 2;
  CSCOM = CSCOM * 2 * chromaRange / colorwheelRes;
  
  if showLabSpace
    figure;
    image(wheelIm);
    blankAxes
    set(gca, 'YDir', 'normal');
    
    hold on;
    plot(imCSCOM(2), imCSCOM(1), 'w+');
    rad = colorwheelRes * chromaRange / overChromaRange / 2; % / colorwheelRes;
    rectangle('Position', [imCSCOM(2)-rad imCSCOM(1)-rad 2*rad 2*rad], 'Curvature', [1 1], 'EdgeColor', 'w');
  end
  
  cform = makecform('lab2srgb');
end


%% Plot points colored by membership

blankFigure;

if cmapByCentroid == 0
  % Use jet
  cmap = colormap(jet(nClusters));
  cmap = [0 0 0; cmap];
else
  % Color the points by the centroids' locations on the usual color wheel
  mu = model.mu;
  
  % Recenter low-D points based on range
  theMin = quantile(pts, colorRangeQuantile);
  theMax = quantile(pts, 1-colorRangeQuantile);
  center = mean([theMin; theMax]);
  mu = bsxfun(@minus, mu, center);
  
  % Renormalize so that the colorScaleQuantile'th point has a distance of 1
  radii = sqrt(sum(pts .^ 2, 2));
  mu = mu / quantile(radii, colorScaleQuantile);
  
  if useLABColor == 0
    cmap = projectOntoColorCone(mu);
  else
    cmap = projectFromLab(mu, LLevel, cform, CSCOM, chromaRange);
  end

  cmap = [0.5 0.5 0.5; cmap];
  colors = cmap(clusters+1, :);
  
  cmap(1, :) = [0 0 0];
end

for c = 0:nClusters
  members = (clusters == c);
  plot(pts(members, 1), pts(members, 2), 'o', 'MarkerFaceColor', cmap(c+1, :), ...
    'MarkerSize', 6, 'MarkerEdgeColor', 'none');
end

axis tight;
axis square;






function wheelImLab = createLABIm(LLevel, chromaRange, colorwheelRes)

LStar = repmat(LLevel, [colorwheelRes colorwheelRes]);
[aStar, bStar] = meshgrid(linspace(-chromaRange, chromaRange, colorwheelRes));
wheelImLab = cat(3, LStar, aStar, bStar);



function map = projectFromLab(pts, LLevel, cform, CSCOM, chromaRange)

pts = fliplr(pts);

% Pull in points that are too far away
angles = atan2(pts(:, 1), pts(:, 2));
radii = sqrt(sum(pts .^ 2, 2));
radii(radii > 1) = 1;
pts(:, 1) = radii .* cos(angles);
pts(:, 2) = radii .* sin(angles);

pts = pts * chromaRange;
pts = bsxfun(@plus, pts, CSCOM);

map = NaN(size(pts, 1), 3);
NaNs = isnan(pts(:, 1));
pts = pts(~NaNs, :);
triples = [repmat(LLevel, size(pts, 1), 1) pts];
map(~NaNs, :) = applycform(triples, cform);


