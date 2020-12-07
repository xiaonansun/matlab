function [movieMC, badFrames, pmtOffFrames] = preprocessCaMovies(tifList, regTif, regFrameNums, dftRegChannel, channels2write, outSuffix, maxMaskWidth, analysisDir, pmt_th, makeMCMrepMovie, frsExclude)
% [movieMC, badFrames, outputsDFT] = preprocessCaMovies_fn2...
%   (tifList, regTif, regFrameNums, dftRegChannel, channels2write [, outSuffix] [, maxMaskWidth] [, analysisDir])
% 
% Takes the names of tif files and information about how to motion correct,
% then performs motion correction, masks off the edges of the movie (which
% are contaminated by motion correction), and returns the movie and which
% frames required too much motion correction.
%
% The motion correction algorithm used is the upsampled discrete Fourier
% transform algorithm in dftregistration.m
%
% INPUTS
%   tifList      -- cell array of strings containing the full path and
%                   filename of each .tif to include
%   regTif       -- string containing the .tif to use for the registration
%                   image (for motion correction)
%   regFrameNums -- vector of frame numbers to use from regTif in making
%                   the registration image
%   dftRegChannel -- channel on which dft registration (motion correction
%                   algorithm) is performed. (same dftoutputs will be used
%                   for registering otherChannels.)
%   channels2write -- optional. Movies of these channels will be written to
%                   tif files. If not provided, all saved channels will be
%                   written to tif.
%   outSuffix    -- optional. If present and non-empty, the resulting movie
%                   will be saved to a series of .tif files, corresponding
%                   to the original .tif files, but with this suffix added
%                   to the filenames. Recommended: '_MCM' (for "motion
%                   corrected and masked"). Default: empty (no saving)
%   maxMaskWidth -- optional. The mask will not be wider than this number
%                   of pixels.
%   analysisDir  -- optional, directory for saving analysis files.
%                   default: same directory as the imaging folder.
%
%
% OUTPUTS
%   movieMC      -- Cell array, each element corresponds to a channel, and
%                   contains the motion corrected and masked movie for that
%                   channel, including frames from all .tif files given.
%   badFrames    -- Cell array, each element corresponds to a channel.
%                   badFrames{ch} is a logical vector of length nFrames,
%                   indicating whether each frame was motion corrected by
%                   more than maxMaskWidth. nFrames is the total number of
%                   frames of all tif files given.
%   pmtOffFrames -- Cell array, each element corresponds to a channel. 
%                   pmtOffFrames{ch} is a logical array indicating frames
%                   during which PMT was off.

%% Optional arguments

if ~exist('outSuffix', 'var')
  outSuffix = '';
end

if ~exist('maxMaskWidth', 'var')
  maxMaskWidth = 20;
end

if ~exist('analysisDir', 'var')
    fPath = fileparts(tifList{1});
else
    fPath = analysisDir;
end

if ~exist('makeMCMrepMovie', 'var')
    makeMCMrepMovie  = 0;
end

if ~exist('frsExclude', 'var')
    frsExclude = [];
end


%% Set mat file name that includes outputsDFT, badFrames, etc.

% Codes below allow for names like "151102_001_002" in case multiple sessions (ie mdf files, aka tif majors) are being analyzed at the same time.
[~, f] = cellfun(@fileparts, tifList, 'uniformoutput', 0);
c = cellfun(@(x)simpleTokenize(x,'_'), f, 'uniformoutput',0);
u = str2num(cell2mat(unique(cellfun(@(x)x{2}, c', 'uniformoutput',0))));
r = repmat('%03d-', 1, length(u)); 
r(end) = [];
d = str2num(cell2mat(unique(cellfun(@(x)x{1}, c, 'uniformoutput',0))));
date_major = sprintf(['%06d_', r], d, u'); % Allows for names like "151102_001_002" in case multiple sessions (ie mdf files, aka tif majors) are being analyzed at the same time.

fileName = fullfile(fPath, [date_major, '.mat']);


%% Get the registration image 
regImage = makeCaImagingRegImage(regTif, regFrameNums, dftRegChannel);

imWidth = size(regImage{dftRegChannel}, 2);
imHeight = size(regImage{dftRegChannel}, 1);


%% Figure out number of recorded channels and how big the resulting movie will be.

nFramesPerMovie = NaN(1, length(tifList));
tifInfo = cell(1, length(tifList));

for t = 1:length(tifList)
    fprintf('Reading info of tif file: %s\n', tifList{t});
    tifInfo{t} = imfinfo(tifList{t});
    channelsSaved = setImagedChannels(tifInfo{t}(1));

    nFramesPerMovie(t) = length(tifInfo{t}) / length(channelsSaved);
end

totalFrames = sum(nFramesPerMovie);


%% To assess motion correction make raw and motion-corrected movies by picking 100 random frames that span the entire length of the movie.

if makeMCMrepMovie     
    nMovFrs = 100; % number of frames of the representative movie.    
    
    %%% Pick random frames that represent the entire movie and append them to imfilename
    ints = floor(totalFrames/nMovFrs);
    rg = 0:ints:totalFrames;
    rg(end) = totalFrames;

    randFrs_motCorrRep = nan(1, nMovFrs);
    for i = 1:nMovFrs
        randFrs_motCorrRep(i) = rg(i)+1 + randi(ints);
    end
    cs0 = [0 cumsum(nFramesPerMovie)] + 1;
    [nrfrs, ~, bn] = histcounts(randFrs_motCorrRep, cs0);
    cs_nrfrs = [0 cumsum(nrfrs)];

    % Initialize representative movies
    movieRawRep = cell(1,max(channelsSaved));
    movieMCMRep = cell(1,max(channelsSaved));
    for ch = channelsSaved
        movieRawRep{ch} = uint16(zeros(imHeight, imWidth, nMovFrs));
        movieMCMRep{ch} = uint16(zeros(imHeight, imWidth, nMovFrs));
    end
    
else
    randFrs_motCorrRep = [];
end


%% Motion correct each movie, gather them together

movieMC = cell(1,max(channelsSaved));
outputsDFT = cell(1,max(channelsSaved));
pmtOffFrames = cell(1,max(channelsSaved)); % frames during which PMT was off. Their pixel values will be set to NaN in movieMC.
for ch = channelsSaved
  movieMC{ch} = uint16(zeros(imHeight, imWidth, totalFrames));
  outputsDFT{ch} = NaN(totalFrames, 4);
  pmtOffFrames{ch} = false(totalFrames, 1);
end
otherChannels = channelsSaved(~ismember(channelsSaved, dftRegChannel)); % dftregistration is done only for dftRegChannel. The same dftoutputs are used for registering otherChannels.

frame = 0;
for t = 1:length(tifList)

    %     fprintf('Motion correcting file: %s\n\n', tifList{t});
    frames = frame + 1 : frame + nFramesPerMovie(t);
    
    frsRep = [];
    if makeMCMrepMovie
        frsRep = randFrs_motCorrRep(bn==t) - cs0(t) + 1; % frame numbers in randFrs_motCorrRep that belong to tif file t
        disp(frsRep)
    end

    
    %% Apply dft registration on dftRegChannel.  
    
    if ~iscell(pmt_th) % a threshold is provided; if frame average intensity is below it, frame will be set as pmtOffFrame.
        [regMovie, DFToutputs, movRep1, pmtOffFrs] = motionCorrectCaImagingFile(tifList{t}, regImage, tifInfo{t}, dftRegChannel, pmt_th, [], 1, frsRep, date_major); 
    else % pmt_th is a cell which contains the index of pmtOffFrames.
        [regMovie, DFToutputs, movRep1] = motionCorrectCaImagingFile(tifList{t}, regImage, tifInfo{t}, dftRegChannel, pmt_th, [], 1, frsRep, date_major); 
        pmtOffFrs = pmt_th{1};
    end
    
    
    %% Save DFToutputs for each tif file.   
    
    [~, fStem] = fileparts(tifList{t});    
    save(fullfile(fPath, fStem), 'DFToutputs')
    
    
    %% Use dftoutputs computed on dftRegChannel for registering otherChannels, if they exist.
    
    if ~isempty(otherChannels)
        [regMovieOther, movRep2] = motionCorrectCaImagingFileCh2(DFToutputs, dftRegChannel, tifList{t}, tifInfo{t}, otherChannels, pmt_th, 1, frsRep, date_major);
        % below will compute different pmtOffFrames for channel2, so we don't use it!
%         [regMovieOther, pmtOffFrsOther] = motionCorrectCaImagingFileCh2(DFToutputs, dftRegChannel, tifList{t}, tifInfo{t}, otherChannels, pmt_th);
    end

    
    %% Gather frames from all tif files.
    
    for ch = dftRegChannel
        movieMC{ch}(:, :, frames) = regMovie{ch};
        outputsDFT{ch}(frames, :) = DFToutputs{ch};
        pmtOffFrames{ch}(frames, :) = pmtOffFrs{ch};
%         pixelShifts{ch}(frames, :) = DFToutputs{ch}(:, 3:4);
    end        
    
    if ~isempty(otherChannels)
        for ch = otherChannels
            movieMC{ch}(:, :, frames) = regMovieOther{ch};
            outputsDFT{ch}(frames, :) = DFToutputs{dftRegChannel};
            pmtOffFrames{ch}(frames, :) = pmtOffFrs{dftRegChannel}; %pmtOffFrsOther{ch}; % Both channels are recorded simultaneously so pmtOffFrames of channel 1 (dftRegChannel) will be used for channel 2 (otherChannels) as well.        
        end
    end
    
    clear regMovie DFToutputs regMovieOther
    
    frame = frame + nFramesPerMovie(t);

    
    %% Save the representative raw movie to compare with its motion-corrected version
    
    if makeMCMrepMovie && ~isempty(cs_nrfrs(t)+1:cs_nrfrs(t+1))
        for ch = dftRegChannel
            movieRawRep{ch}(:, :, cs_nrfrs(t)+1:cs_nrfrs(t+1)) = movRep1{ch};
        end
        
        if ~isempty(otherChannels)
            for ch = otherChannels
                movieRawRep{ch}(:, :, cs_nrfrs(t)+1:cs_nrfrs(t+1)) = movRep2{ch};
            end
        end
    end
    
end


%% Mask result, to get rid of edges that are affected by motion correction

maskBounds = cell(1,max(channelsSaved));
badFrames = cell(1,max(channelsSaved));
for ch = channelsSaved
    if maxMaskWidth > 0
        
      % maskBounds shows the bounds (left, right, top and bottom) of motion-corrected movie relative to
      % the original, border-trimmed movie (size 512x402).
      [maskBounds{ch}, badFrames{ch}] = determineMovieMaskBounds(outputsDFT{ch}(:, 3:4), [imWidth imHeight], maxMaskWidth, 0, frsExclude); % maskBounds shows bad pixels on left, right, top and bottom, respectively; so it corresponds to columns and rows, respectively.

      fprintf('Masking off pixels (channel %d): %d on left, %d on right, %d on top, %d on bottom\n', ...
        ch, maskBounds{ch}(1)-1, imWidth-maskBounds{ch}(2), maskBounds{ch}(3)-1, imHeight-maskBounds{ch}(4));

    
        movieMC{ch} = maskMovie(movieMC{ch}, maskBounds{ch});

    else        
      badFrames{ch} = true(size(movieMC{ch}, 3), 1);
    end
end


%% save outputsDFT, badFrames and maskBounds corresponding to all frames (ie all tif parts) of a mdf file.

%{
[~, fStem] = fileparts(tifList{1});
tokens = simpleTokenize(fStem,'_'); 
fileName = fullfile(fPath, [tokens{1}, '_', tokens{2}, '.mat']);
%}

if exist(fileName, 'file')==2 % date_major mate file alrady exists, so append to it!
    save(fileName, '-append', 'outputsDFT', 'badFrames', 'pmtOffFrames', 'maskBounds')
else
    save(fileName, 'outputsDFT', 'badFrames', 'pmtOffFrames', 'maskBounds')
end


%% save badFrames and pmtOffFrames for each tif movie separately.

cs = [0 cumsum(nFramesPerMovie)];
for t = 1:length(tifList)
    
    [~, fStem] = fileparts(tifList{t});  
    frames = cs(t)+1 : cs(t+1);
    badFramesTif = cell(1,max(channelsSaved));
    pmtOffFramesTif = cell(1,max(channelsSaved));    
    
    for ch = 1:max(channelsSaved)
        if ~isempty(badFrames{ch})
            badFramesTif{ch} = badFrames{ch}(frames);
            pmtOffFramesTif{ch} = pmtOffFrames{ch}(frames);
        end
    end
    save(fullfile(fPath, fStem), 'badFramesTif', 'pmtOffFramesTif', '-append')
end
    
clear badFramesTif pmtOffFramesTif


%% Write movie to file, if requested

if ~isempty(outSuffix)
    
    if ~exist('channels2write', 'var')
        channels2write = channelsSaved;
    end
    
    for ch = channels2write % unlike the raw tifs, the motion corrected tifs will have channel 1 and channel 2 in separate files, instead of alternating frames in the same file.
      % The apparent brightness changes, but I think this is just a scaling issue
      % from a header parameter I can't change

      frame = 0;
      for t = 1:length(tifList)
          
        % Figure out filename
        [~, fStem, fExt] = fileparts(tifList{t});
%         outFile = fullfile(fPath, [fStem outSuffix fExt]);
        outFile = fullfile(fPath, [[fStem,'_ch',num2str(ch)] outSuffix fExt]);

        fprintf('Writing file %s (%d/%d)\n', outFile, t, length(tifList));

        % Figure out frames
        frames = frame + 1 : frame + nFramesPerMovie(t);

        imwrite(movieMC{ch}(:, :, frames(1)), outFile, 'TIF', ...
          'Resolution', [size(movieMC{ch}, 2) size(movieMC{ch}, 1)], 'Compression', 'none');

        if length(frames) > 1
          for f = 2:length(frames)
            if mod(f, 100) == 0
              fprintf('%d ', f);
            end
            if mod(f, 1000) == 0
              fprintf('\n');
            end

            imwrite(movieMC{ch}(:, :, frames(f)), outFile, 'TIF', ...
              'Resolution', [size(movieMC{ch}, 2) size(movieMC{ch}, 1)], 'Compression', 'none', ...
              'WriteMode', 'append');
          end
        end

        % Save this chunk of badFrames and pixelShifts
%         badFramesTif = badFrames{ch}(frames);
%         pixelShiftsTif = pixelShifts(frames, :);
%         save(fullfile(fPath, [fStem '_badFrames']), 'badFramesTif', 'pixelShiftsTif');    
    
        frame = frame + nFramesPerMovie(t);
      end
      
    end
end

fprintf('\n---------- Done ----------.\n\n\n');
    

%% Ser the representative motion-corrected movie and save the vars

if makeMCMrepMovie   
%     frames = randFrs_motCorrRep;
    for ch = 1:length(movieMC)
        if ~isempty(movieMC{ch})
            movieMCMRep{ch} = movieMC{ch}(:, :, randFrs_motCorrRep);
            %{
            temp = movieMC{ch}(:, :, frames);
            outFile = sprintf('movRepMCM_ch%i', ch);
            eval([outFile '= temp;'])
            save(fileName, outFile, '-append')
            %}
        end
    end
    save(fileName, 'randFrs_motCorrRep', 'movieMCMRep', 'movieRawRep', '-append')
end


    