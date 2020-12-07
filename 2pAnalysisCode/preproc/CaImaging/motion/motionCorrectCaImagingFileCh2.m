function [regMovie, movRep, pmtOffFrsOther] = motionCorrectCaImagingFileCh2(outputsDFT, dftRegChannel, tifName, tifInfo, otherChannels, pmt_th, trimBorders, randFrs_motCorrRep, date_major)
%
% This funciton uses dftoutputs computed on dftRegChannel to regist
% otherChannels.
%
% regMovie  -- Cell array, each element contains the motion corrected 
%              movie of a channel defined in otherChannels.
% pmtOffFrsOther -- Cell array, each element corresponds to a channel. pmtOffFrsOther{ch} is a
%                   logical array indicating frames during which PMT was off.

%% Parameters

% Optional arguments
if ~exist('trimBorders', 'var')
    trimBorders = 1;
end

if ~exist('randFrs_motCorrRep', 'var')
    randFrs_motCorrRep = [];
end


% This is how wide the black borders on the left and right sides of the
% image are, when using the MScan option to correct for the sinusoidal
% movement of the mirrors and a horizontal resolution of 512. These borders
% will get chopped off.
borderWidth = 55;

channelsSaved = setImagedChannels(tifInfo(1));

nFrames = length(tifInfo);
imWidth = tifInfo(1).Width;
imHeight = tifInfo(1).Height;

% Prepare to trim borders
if trimBorders
    validPixels = [false(1, borderWidth) true(1, imWidth - 2*borderWidth) false(1, borderWidth)];
else
    validPixels = true(1, imWidth);
end


%%%
% if ~isempty(randFrs_motCorrRep)
movRep = cell(1,max(channelsSaved));
% end


%%
regMovie = cell(1,max(channelsSaved));
pmtOffFrsOther = cell(1,max(channelsSaved));

tic;
for ch = otherChannels;
    
    framesToUse = ch : length(channelsSaved) : nFrames;
    
    %% Read frames corresponding to otherChannels out of the tiff and trim borders.
    
    fprintf('Reading tiff %s, channel %d\n', tifName, ch);
    
    movie = bigread2_frs2use(tifName, framesToUse);
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
    
    if nargout > 2
        pmtOffFrsOther{ch} = false(1,size(movie,3));
        if ~isnan(pmt_th)
            framesAve = mean(reshape(movie, size(movie,1)*size(movie,2), []));        
            pmtOffFrsOther{ch}(framesAve < pmt_th) = true;

            movie(:,:,pmtOffFrsOther{ch}) = NaN;
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
    end
    

    %% Motion correction / registration
    
    fprintf('Correcting motion, channel %d\n', ch);
    
    for f = 1:size(movie, 3)
        % Display progress
        if mod(f, 100) == 0
            fprintf('%d ', f);
        end
        if mod(f, 1000) == 0
            fprintf('\n');
        end

        %% use the same outputsDFT computed on dftRegChannel for otherChannels
        buf2ft = fft2(movie(:, :, f));
        %     output=[error,diffphase,row_shift,col_shift];
        diffphase = outputsDFT{dftRegChannel}(f,2);
        row_shift = outputsDFT{dftRegChannel}(f,3);
        col_shift = outputsDFT{dftRegChannel}(f,4);
        
        
        %% register movie
        [nr,nc] = size(buf2ft);
        Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
        Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
        [Nc,Nr] = meshgrid(Nc,Nr);
        Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
        Greg = Greg*exp(i*diffphase);
        
        
        %%
        regMovie{ch}(:, :, f) = uint16(abs(ifft2(Greg)));
        
    end
    fprintf('\n');
    
end

fprintf('Registering %d frames from %d channel(s) took %0.1f s\n\n', length(framesToUse), length(otherChannels), toc);


