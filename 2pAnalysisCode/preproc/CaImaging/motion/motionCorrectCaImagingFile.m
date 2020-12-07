function [regMovie, outputsDFT, movRep, pmtOffFrs] = motionCorrectCaImagingFile(tifName, refImage, tifInfo, dftRegChannel, pmt_th, framesToUse, trimBorders, randFrs_motCorrRep, date_major)
% [regMovie, outputsDFT] = motionCorrectCaImagingFile(tifName, refImage, tifInfo, dftRegChannel [, framesToUse] [, trimBorders])
%
% Use the efficient subpixel registration algorithm to motion correct a
% calcium imaging movie.  Uses 10-fold upsampling.
% Check this link for more info: http://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation/content/html/efficient_subpixel_registration.html
%
% INPUTS
%   tifName      -- name of the file to operate on
%   refImage     -- reference image to register to (get from makeCaImagingRegImage() )
%   tifInfo      -- imfinfo of the file to operate on.
%   dftRegChannel -- numeric array of channel(s) to use for dft registration. 
%                   If only one channel was saved, that channel will be used
%                   for motion correction and dftRegChannel is not needed.
%   framesToUse  -- optional. If provided and not NaN, will extract only
%                   these frames
%   trimBorders  -- optional (default true). If using the MScan option to
%                   correct for sinusoidal motion of the fast mirror, there
%                   will be black bars on the left and right of the image.
%                   This option lets you throw those away (assumes 55 pixel
%                   border width, which is correct for a 512 pixel wide
%                   image).
%
% OUTPUTS
%   regMovie     -- Cell array, each element contains the motion corrected 
%                   movie of a channel defined in dftRegChannel. (Each frame 
%                   of regMovie is the raw image after being shifted to get 
%                   aligned with the reference image.)
%   outputsDFT   -- Cell array, each element corresponds to a channel 
%                   defined in dftRegChannel. outputsDFT{ch} has size
%                   nFrames x 4, and contains the 1st output of
%                   dftregistration = [error, diffphase, net_row_shift,
%                   net_col_shift], . The last 2 columns show pixel shifts
%                   in [y x]. (See dftregistration for more details.)
%   pmtOffFrs    -- Cell array, each element corresponds to a channel. pmtOffFrs{ch} is a
%                   logical array indicating frames during which PMT was off.
%
% To trim off edges after motion correction, use maskMovie()


%% Parameters

% This is how wide the black borders on the left and right sides of the
% image are, when using the MScan option to correct for the sinusoidal
% movement of the mirrors and a horizontal resolution of 512. These borders
% will get chopped off.
% Note: this is now detected automatically, but the below makes it so that
% for images with very dark edges we still get a correct answer. This works
% because the software has always used one or the other of these two
% values.
validBorderWidths = [42 55];

% Upsampling factor for subpixel motion correction. 10 seems likes more
% than enough.
usFac = 10;

framesToUseDefault = false;

% Optional arguments
if ~exist('framesToUse', 'var') || isempty(framesToUse)
  framesToUseDefault = true; % if no input is provided all frames will be analyzed.
end

if ~exist('trimBorders', 'var')
  trimBorders = 1;
end

if ~exist('randFrs_motCorrRep', 'var')
    randFrs_motCorrRep = [];
end
   

%% Read tiff metadata

% tifInfo = imfinfo(tifName);
channelsSaved = setImagedChannels(tifInfo(1));
if length(channelsSaved)==1
    dftRegChannel = channelsSaved;
end


nFrames = length(tifInfo);
imWidth = tifInfo(1).Width;
imHeight = tifInfo(1).Height;

% Prepare to trim borders
validPixels = [];
% if trimBorders
%   validPixels = [false(1, borderWidth) true(1, imWidth - 2*borderWidth) false(1, borderWidth)];
% else
%   validPixels = true(1, imWidth);
% end
if ~trimBorders
  validPixels = true(1, imWidth);
end


%%%
% if ~isempty(randFrs_motCorrRep)
movRep = cell(1,max(channelsSaved));
% end


%%
regMovie = cell(1,max(channelsSaved));
outputsDFT = cell(1,max(channelsSaved));
pmtOffFrs = cell(1,max(channelsSaved));

tic;
for ch = dftRegChannel
    
    if framesToUseDefault
        framesToUse = find(ismember(channelsSaved, dftRegChannel)) : length(channelsSaved) : nFrames;
    end

    %% Read all the images out of the tiff and trim borders.
    
    fprintf('Reading tiff %s, channel %d\n', tifName, ch);
    
    movie = bigread2_frs2use(tifName, framesToUse);
    
    % If not found already, figure out how wide the left and right black
    % borders of the image are
    if trimBorders && isempty(validPixels)
      meanMov = mean(movie(:, :, 1:min(100, size(movie, 3))), 3);
      margMov = mean(meanMov, 1);
      
      % If not found already, figure out how wide the left and right black
      % borders of the image are.
      % I needed to make this robust, because sometimes the registration
      % image will have such a dark edge that we mis-estimate the borders. So
      % now we round off to one of two values, which were what got used in
      % different versions of the MScan software.
      LBWidth = find(margMov > 0, 1) - 1;
      RBWidth = imWidth - find(margMov > 0, 1, 'last');
      bWidthEst = min([LBWidth RBWidth]);
      [minBorderError, probableBorderI] = min(abs(validBorderWidths - bWidthEst));
      bWidth = validBorderWidths(probableBorderI);
      if minBorderError > 0
        warning('Possible border detection error: guess was %d, correcting to %d', bWidthEst, bWidth);
      end
      validPixels = [false(1, bWidth) true(1, imWidth - 2 * bWidth) false(1, bWidth)];
    end
    
    movie = movie(:, validPixels, :);
    
    %{
    % Pre-allocate movie
    movie = zeros(imHeight, sum(validPixels), length(framesToUse), 'uint16');

    % Read frames, throwing away borders
    for f = 1:length(framesToUse)
      if mod(f, 100) == 0
        fprintf('%d ', f);
      end
      if mod(f, 1000) == 0
        fprintf('\n');
      end  
      rawFrame = imread(tifName, 'Index', framesToUse(f), 'Info', tifInfo);
      movie(:, :, f) = rawFrame(:, validPixels);
    end
    %}
    
    %%
    clear rawFrame
    
    %% Find the frames recorded while PMT was off, and set their pixel values to NaN. (Note: since movie is uint16, NaN will be turned to 0)
    
    if nargout > 3
        pmtOffFrs{ch} = false(1,size(movie,3));
        if ~isnan(pmt_th)
            framesAve = mean(reshape(movie, size(movie,1)*size(movie,2), []));        
            pmtOffFrs{ch}(framesAve < pmt_th) = true;

            movie(:,:,pmtOffFrs{ch}) = NaN;
        end
    end
    

    %% Save a movie of random frames to serve as a raw movie example (to compare it with the motion-corrected one if desired).

    if ~isempty(randFrs_motCorrRep)
        movRep{ch} = movie(:, :, randFrs_motCorrRep);
        
        %{
        frames = randFrs_motCorrRep;                    
        temp = movie(:, :, frames);
        outFile = sprintf('movRepRaw_ch%i', ch);
        eval([outFile '= temp;'])            

        save(fullfile(fileparts(tifName), date_major), outFile, '-append')
        %}
        % make a tif file
        %{
        imwrite(movieMC{ch}(:, :, frames(1)), outFile, 'TIF', ...
            'Resolution', [size(movieMC{ch}, 2) size(movieMC{ch}, 1)], 'Compression', 'none');

        if length(frames) > 1
            for f = 2:length(frames)                
                imwrite(movieMC{ch}(:, :, frames(f)), outFile, 'TIF', ...
                    'Resolution', [size(movieMC{ch}, 2) size(movieMC{ch}, 1)], 'Compression', 'none', ...
                    'WriteMode', 'append');
            end
        end        
        %}        
%     else
%       movRep = {};
    end
    
    
    %% Motion correction / registration

    fprintf('Correcting motion, channel %d\n', ch);

    % Get FFT of reference registration image
    fftRef = fft2(refImage{ch});

    % Pre-allocate result
    regMovie{ch} = uint16(zeros(size(movie)));
    outputsDFT{ch} = zeros(size(movie, 3), 4);

    % Pre-allocate for saving the shift magnitudes
    dftOutputs = zeros(size(movie, 3), 4);

    % Do the registration
%     tic;
    for f = 1:size(movie, 3)
      % Display progress
      if mod(f, 100) == 0
        fprintf('%d ', f);
      end
      if mod(f, 1000) == 0
        fprintf('\n');
      end

      [dftOutputs(f, :), Greg] = dftregistration(fftRef, fft2(movie(:, :, f)), usFac);
      regMovie{ch}(:, :, f) = uint16(abs(ifft2(Greg)));
    end
    fprintf('\n');
%     fprintf('\n\nRegistering %d frames took %0.1f s\n\n', size(movie, 3), toc);

    %%% pixelsShifts output
%     pixelShifts{ch} = dftOutputs(:, 3:4);
    outputsDFT{ch} = dftOutputs;

    
    %%
    clear fftRef dftOutputs Greg movie
    
end

fprintf('Registering %d frames from %d channel(s) took %0.1f s\n\n', length(framesToUse), length(dftRegChannel), toc);

    
