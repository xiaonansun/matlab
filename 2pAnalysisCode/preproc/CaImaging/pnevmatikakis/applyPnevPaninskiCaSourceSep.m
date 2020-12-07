function [A, C, b, f, S, Df, C_df, S_df, P, ROIOrderMeasure, Ain] = applyPnevPaninskiCaSourceSep(movie, userParams)
% 
% [spatialBasis, activity, spatialBackground, backgroundActivity, spiking, Df, activityDf, spikingDf, pnevParam, ROIOrderMeasure, greedySpatialBasis] = ...
%    applyPnevPaninskiCaSourceSep(movie [, userParams])
%
% INPUTS
%   movie      -- the height x width x frames movie. If it is not in double
%                 format, it will be automatically converted.
%   userParams -- optional struct. Parameters to affect the algorithm. Each
%                 field is also optional. Possible fields are:
%      .maxGreedy  -- maximum number of frames to be processed for greedy
%                     spatial initialization. Default 2000.
%      .nComps     -- number of ROIs (components) to be found. Default 300.
%      .gSiz       -- maximum size of initial spatial footprint (box of
%                     size gSiz x gSiz). Default 15.
%      .gSig       -- typical size of neuron, in std dev of gaussian.
%                     Default 3.
%      .makePlots  -- whether or not to make some basic plots. Should be 0
%                     if using -nodesktop. This feature is not very useful,
%                     and might be removed in the future. Default 0.
%      .makeMovies -- whether or not to make some basic diagnostic movies.
%                     Should be 0 if using -nodesktop. This feature is not
%                     very useful, and might be removed in the future.
%                     Default 0.
%
% OUTPUTS
%   spatialBasis       -- nPixels x nROIs. Where the neurons are.
%   activity           -- nROIs x nFrames. Inferred activity of each ROI
%                         during each frame.
%   spatialBackground  -- nPixels x 1. Spatial structure of the background
%                         activity (effectively a picture of the neuropil).
%   backgroundActivity -- 1 x nFrames. Timecourse of background (neuropil)
%                         activity.
%   pixelNoise         -- nPixels x 1. How much each pixel varied with time
%   activityDF         -- like activity, but extracted to dF/F using
%                         Eftychios's method. Needs comparison with the
%                         extraction I usually do, using a variant of
%                         Arthur Konnerth's method.
%   spiking            -- nROIs x nFrames inferred spike events. Not sure
%                         what the units of measurement are.
%   ROIOrderMeasure           -- nROIs x 1. How "good" each ROI is likely to be.
%                         Will be sorted to be monotonically decreasing.
%                         From Eftychios's order_ROIs(), using L4 norm of
%                         max of each pixel (I think).
%
% based on Eftychios's demo_script.m from 10/5/15


%% Parameters

% For greedyROI
params.maxGreedy = 2000;
params.nComps = 300;
params.gSiz = 15;
params.gSig = 3;
params.makePlots = 0;
params.makeMovies = 0;


%% Override parameters with user-provided parameters if supplied

if exist('userParams', 'var')
  fields = fieldnames(userParams);
  
  for fi = 1:length(fields)
    if isfield(params, fields{fi})
      params.(fields{fi}) = userParams.(fields{fi});
    else
      fprintf('applyPnevPaninskiCaSourceSep: UNKNOWN PARAMETER: %s\n', fields{fi});
    end
  end
end


%% Convert movie data type to double, prep sizing info

if ~isa(movie, 'double');
  fprintf('Converting movie format to double... ');
  tic;
  movie = double(movie);
  fprintf('%0.2fs\n', toc);
end

[imHeight, imWidth, nFrames] = size(movie);
nPixels = imHeight * imWidth;


%% Could insert a missing data interpolation step here, if needed


%% Fast initialization of spatial components using the greedyROI
% greedyROI2d takes a greedy guess at where ROIs are
% This is using a subsample of the movie, to save memory and perform
% faster.

fprintf('Initializing with greedy algorithm...\n');
tic;

% Ain is the spatial components, ie the set of ROI masks (one per neuron),
% Cin is the temporal components, ie the initial guess at the fluorescence
% trace for each neuron based on these masks, and center is the center of
% each ROI
if nFrames > params.maxGreedy
  frames = round(linspace(1, nFrames, params.maxGreedy));
  [Ain, Cin, center, ~] = greedyROI2d(movie(:,:,frames), params.nComps, params);
else
  [Ain, Cin, center, ~] = greedyROI2d(movie, params.nComps, params);
end

% Ain is the initial pixel-vector version of basis; that is, where each ROI is
Ain = sparse(reshape(Ain, nPixels, params.nComps)); % initial estimate of spatial footprints (size d x nr)
Cin = Cin'; % initial estimate of temporal components. (size nr x T)

fprintf('%0.2fs\n', toc);
    

%% Display centers of found components with greedy algorithm

if params.makePlots
    Cn =  correlation_image(movie); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
    figure;
    subplot(121), imagesc(Cn);
    axis equal; axis tight; hold all;
    
    subplot(122), imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;
end


%% Compute estimates of noise for every pixel and a global time constant

% Yr is the reshaped observations (movie)
% The sequence below is a little weird, because I want to rename the
% variable but want Matlab to use in-place operations
movie = reshape(movie, nPixels, nFrames);
Yr = movie;
clear movie;

p = 2;                                                          % order of autoregressive system (p=1 just decay, p=2, both rise and decay)
active_pixels = find(sum(Ain, 2));                              % pixels where the greedy method found activity
unsaturated_pixels = find_unsaturatedPixels(Yr);                % pixels that do not exhibit saturation
options.pixels = intersect(active_pixels, unsaturated_pixels);  % base estimates only on unsaturated

P = arpfit(Yr, p, options);  % estimate noise for each pixel & global time constant


% Perform non-negative matrix factorization of the background activity --
% the movie residuals, which are the movie minus the activity attributed to
% the known neurons. This step factors the background into a spatial matrix
% (bin) and a timecourse matrix (fin), respectively.
if nFrames > params.maxGreedy
  [bin, fin] = nnmf(max(Yr(:, frames) - Ain * Cin, 0), 1);
else
  [bin, fin] = nnmf(max(Yr - Ain * Cin, 0), 1);
end

% Not using missing data interpolation, so not using the below...
% P.interp = Y_interp;
% % remove interpolated values
% miss_data_int = find(Y_interp);
% Yr(miss_data_int) = P.interp(miss_data_int);

P.unsaturatedPix = unsaturated_pixels;


%% If conserving memory in greedy step, estimate full activity time course

if nFrames > params.maxGreedy
  fprintf('Initial temporal component estimation...\n');
  tic;
  
  % Note: at Eftychios's suggestion, using 4 temporal iterations instead of
  % 2 because we don't have a good initial estimate for C
  P.method = 'constrained_foopsi';            % choice of method for deconvolution
  P.temporal_iter = 4;                        % number of iterations for block coordinate descent
  P.fudge_factor = 0.98;                      % fudge factor to reduce time constant estimation bias
  
  % C is the new estimate of temporal components, f is the
  % background timecourse, Yres is the new set of movie residuals, P is the
  % params struct (which is augmented inside the function)
  [Cin, fin, ~, P, ~, Ain] = update_temporal_components(Yr, Ain, bin, [], [], P);

  fprintf('%0.2fs\n', toc);
end


%% Update spatial components

fprintf('Updating spatial components...\n');
tic;

P.search_method = 'ellipse';
P.d1 = imHeight;
P.d2 = imWidth;
P.dist = 3;                 % ellipse expansion factor for local search of spatial components

% A is the updated spatial layouts of the ROIs, b is the updated spatial
% distribution of the background
[A, b] = update_spatial_components(Yr, Cin, fin, Ain, P);
% [A, b, C] = update_spatial_components(Yr, Cin, fin, Ain, P);

fprintf('%0.2fs\n', toc);


%% Update temporal components

fprintf('Updating temporal components...\n');
tic;

P.method = 'constrained_foopsi';            % choice of method for deconvolution
P.temporal_iter = 2;                        % number of iterations for block coordinate descent
P.fudge_factor = 0.98;                      % fudge factor to reduce time constant estimation bias

% C is the new estimate of temporal components, f is the updated
% background timecourse, Yres is the new set of movie residuals (= Y - A*C
% - b*f; d X T matrix), P is the params struct (which is augmented inside
% the function) 
[C, f, Yres, P, S, A] = update_temporal_components(Yr, A, b, Cin, fin, P);
% [C, f, Yres, P, S, A] = update_temporal_components(Yr, A, b, C, fin, P);

fprintf('%0.2fs\n', toc);


%% Merge found ROIs based on overlap and correlation

% P.method = 'project';
P.merge_thr = 0.9;                           % merging threshold

fprintf('Merging ROIs...\n');

[Am, Cm, nComps, merged_ROIs, P, Sm] = merge_ROIs(Yres, A, b, C, f, P, S);
% [A, C, nComps, mergedROIs, P, S] = merge_ROIs(Yres, A, b, C, f, P, S);

fprintf('Reduced to %d ROIs\n', nComps);


%% Display merging

% display_merging = 0; % flag for displaying merging example
if params.makePlots % display_merging
    for i = 1:length(merged_ROIs) % 1; randi(length(merged_ROIs));
        ln = length(merged_ROIs{i});
        figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),imHeight,imWidth));
            title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,nComps-length(merged_ROIs)+i),imHeight,imWidth));
        title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight;
        
        subplot(1,ln+2,ln+2);
        plot(1:nFrames,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))');
        % plot merged component in black
        hold all; plot(1:nFrames,Cm(nComps-length(merged_ROIs)+i,:)/max(Cm(nComps-length(merged_ROIs)+i,:)),'--k')
        title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
    end
end


%% Repeat

fprintf('Updating temporal components with new ROIs...\n');
tic;

[A2, b2] = update_spatial_components(Yr, Cm, f, Am, P);
[C2, f2, Yres, P, S2, A2] = update_temporal_components(Yr, A2, b2, Cm, f, P);
nComps = size(C2,1); % some components might be removed due to complex values, so reset nComps.
[C_df, Df, S_df] = extract_DF_F(Yr, [A2, b2], [C2;f2], S2, nComps + 1); % extract DF/F values (optional)


% Get rid of nuisance background component that comes out of extract_DF
C_df(end, :) = [];
Df(end) = [];

fprintf('%0.2fs\n', toc);


%% Reorder the ROIs based on their maximum temporal activation and their size (through their l_inf norm)

[A, C, S, P, ~, ROIOrderMeasure] = order_ROIs(A2, C2, S2, P);


%% Final plotting, if desired

if params.makePlots
  contour_threshold = 0.95;             % amount of energy used for each component to construct contour plot
  figure;
  [Coor, json_file] = plot_contours(A, reshape(P.sn, imHeight, imWidth), contour_threshold, 1); % contour plot of spatial footprints
  %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
  
  % The below makes a plot of activity in each patch
  view_patches(Yr, A, C, b2, f2, imHeight, imWidth)
end

if params.makeMovies  
  % Make an .avi movie of activity in a few patches
  param.skip_frame = 2;
  param.ind = [1, 2, 3, 4];
  param.sx = 16;
  param.make_avi = 0;
  make_patch_video(A, C, b2, f2, Yr, imHeight, imWidth, param)
end


%% Prepare remaining outputs

% A = reshape(A, [imHeight imWidth size(A, 2)]);

% Pixel noise
% sn = reshape(P.sn, imHeight, imWidth);
% sn = P.sn;
% g = P.g;

