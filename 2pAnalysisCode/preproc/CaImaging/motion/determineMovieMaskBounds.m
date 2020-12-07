function [bounds, badFrames] = determineMovieMaskBounds(pixelShifts, movieRes, maxMaskWidth, showShiftHisto, frsExclude)
% [bounds, badFrames] = determineMovieMaskBounds(pixelShifts, movieRes, maxMaskWidth [, showShiftHisto])
%
% Choose the bounds for masking the motion-corrected movie. This function
% will only produce movies with height and width being even numbers,
% because Eftychios's algorithm won't extract one row/column of pixels if
% that dimension is odd.
%
% INPUTS
%   pixelShifts   -- output from motionCorrectCaImagingFile()
%   movieRes      -- [width height] of images
%   maxMaskWidth  -- the widest mask permitted. If there are shifts larger
%                    than this, those frames will be marked bad; Set to inf
%                    to allow masking as much as needed.
%   showShiftHisto -- optional. If true, display a histogram of x shifts
%                     and a histogram of y shifts. Default false
%
% OUTPUTS
%   bounds        -- the mask bounds, to use as input to maskMovie. bounds
%                    is [x1 x2 y1 y2]
%   badFrames     -- a logical vector with size nFrames x 1. When a frame
%                    has a larger pixel shift during motion correction than
%                    the width of the mask, it is marked as bad (true).


%% Optional arguments

if ~exist('showShiftHisto', 'var')
  showShiftHisto = 0;
end


%% Find largest pixel shifts in each direction

xShifts = pixelShifts(:, 2);
yShifts = pixelShifts(:, 1);

% If asked, plot histograms of the pixel shifts
if showShiftHisto
  figure;
  hist(xShifts, 50);
  title('x shifts');
  
  figure;
  hist(yShifts, 50);
  title('y shifts');
end


if exist('frsExclude', 'var') % frames not to use for computing the mask; they will be set to badFrames.                                         
    xShifts(frsExclude) = nan;
    yShifts(frsExclude) = nan;
end

% How much images were shifted in each direction (rounded up)
xShiftMin = floor(min(xShifts));
xShiftMax = ceil(max(xShifts));
yShiftMin = floor(min(yShifts));
yShiftMax = ceil(max(yShifts));


%% Apply maxMaskWidth limit to shifts
% That is, never mask off more then maxMaskWidth pixels on a side

xShiftMin = max([-maxMaskWidth xShiftMin]);
xShiftMax = min([maxMaskWidth xShiftMax]);
yShiftMin = max([-maxMaskWidth yShiftMin]);
yShiftMax = min([maxMaskWidth yShiftMax]);


%% Ensure shifts are in the right direction
% This will generally be the case already, but might not be if the
% reference image is from a different file

xShiftMin = min([xShiftMin 0]);
xShiftMax = max([xShiftMax 0]);
yShiftMin = min([yShiftMin 0]);
yShiftMax = max([yShiftMax 0]);


%% Make sure bounds give us width and height that are multiples of 2

% Clip extra pixel for width if needed
if mod(movieRes(1) + xShiftMin - xShiftMax, 2) ~= 0
  % Need extra logic so we don't trim back more than the maxMaskWidth
  if -xShiftMin < maxMaskWidth
    xShiftMin = xShiftMin - 1;
  else
    xShiftMax = xShiftMax + 1;
  end
end

% Clip extra pixel for height if needed
if mod(movieRes(2) + yShiftMin - yShiftMax, 2) ~= 0
  % Need extra logic so we don't trim back more than the maxMaskWidth
  if -yShiftMin < maxMaskWidth
    yShiftMin = yShiftMin - 1;
  else
    yShiftMax = yShiftMax + 1;
  end
end


%% Determine mask

% A positive shift for x indicates the image had to be shifted right to get
% aligned (ie mouse brain moved left), so clip the left (low index) edge (since
% the algorithm we use for registration uses wraparound). 
% A positive shift for y indicates the image had to be shifted down to
% get aligned (ie mouse brain moved up), so clip the top (low index) edge.
%
% FN: The variable "bounds" indicates [left, right, top, down] bounds on the
% registered movie (512x402 pixels). "bounds" will be used to clip the movie
% edges. The edges need to be removed because they include wrapped-around
% pixels as a result of dft registration.
%
% For a 512x512 image, bounds = [left, right, top, down] corresponds to
% following pixel numbers:
% top minus 1 : # pixels to be clipped on the top 
% 512 minus down : # pixels to be clipped on the bottom 
% left minus 1 : # pixels to be clipped on the left 
% 402 minus right : # pixels to be clipped on the right (Note: it is 402 because we already removed the 55-pixel dark borders on the left and right: 512 - 55*2 = 402).

bounds = [xShiftMax+1 movieRes(1)+xShiftMin yShiftMax+1 movieRes(2)+yShiftMin];


%% Figure out which frames had shifts exceeding mask bounds

% Remember if xShifts (or yShifts) for a frame is more than maxMaskWidth,
% that frame will include some columns (or rows) on the edges that do not
% have any correspondant part in the reference image (since you are not
% clipping more than maxMaskWidth.) 

% This is easy, because if there were any too-big shifts, then the mask
% will be maximum width in that direction

biggerShift = max(abs(pixelShifts), [], 2);
badFrames = (biggerShift > maxMaskWidth);
if exist('frsExclude', 'var')
    badFrames(frsExclude) = true;
end

fprintf('Fraction of badFrames (ie frames with more pixelShifts than maxMaskWidth) = %.2f\n', nanmean(badFrames))

