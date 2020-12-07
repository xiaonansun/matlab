function regImage = makeCaImagingRegImage(tifName, regFrameNums, dftRegChannel, trimBorders)
% regImage = makeCaImagingRegImage(tifName, frames, dftRegChannel [,trimBorders])
%
% Make a median image to use as a reference (registration image) for motion
% correction. tifName should be the name of the tiff file to use, and
% frames should contain a vector of frame numbers (within that file).
%
% If using the MScan option to correct for sinusoidal motion of the fast
% mirror, there will be black bars on the left and right of the image. The
% optional trimBorders argument lets you throw those away (assumes 55 pixel
% border width, which is correct for a 512 pixel wide image). Must agree
% with option used in motionCorrectCaImagingFile(). Default true.

% dftRegChannel: numeric array of channel(s) to apply dft registration on. 
%     If only one channel was saved, that channel will be used 
%     for motion correction and dftRegChannel is not needed.
%
% regImage: cell array, each element includes the regImage of the
%     corresponding channel. e.g. if dftRegChannel = [2], 
%     regImage{1} will be empty, and regImage{2} will be the regImage 
%     of channel 2.


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


%% Optional arguments

if ~exist('trimBorders', 'var')
  trimBorders = 1;
end


%% Read metadata

% fprintf('Reading tiff for registration image\n\n');

tifInfo = imfinfo(tifName);
imWidth = tifInfo(1).Width;
imHeight = tifInfo(1).Height;


%% Channels to use for motion correction

channelsSaved = setImagedChannels(tifInfo(1));

if length(channelsSaved)==1
    dftRegChannel = channelsSaved;
end


%%
% regFrameNums(regFrameNums > length(tifInfo)/length(channelsSaved)) = []; % regTif must be the 1st tif file.
regImage = cell(1, max(channelsSaved));

cnt = 0;
for ch = dftRegChannel;
    cnt = cnt+1;
    frames = regFrameNums*length(channelsSaved) - (length(channelsSaved) - cnt); % frame numbers corresponding to regFrameNums in the tif file (which includes frames from all channels).
  
    %% Read the chosen frames out of the tiff

    % Pre-allocate movie
    movie = zeros(imHeight, imWidth, length(frames), 'uint16');

    fprintf('Reading reference tif file %s, %d frames....\n', tifName, length(frames))
    % Read frames
    for f = 1:length(frames)
      movie(:, :, f) = imread(tifName, 'Index', frames(f), 'Info', tifInfo);
    end
   

    %% Take the median to get regImage

    regImage{ch} = median(movie, 3);


    %% Trim borders

    % If not found already, figure out how wide the left and right black
    % borders of the image are.
    % I needed to make this robust, because sometimes the registration
    % image will have such a dark edge that we mis-estimate the borders. So
    % now we round off to one of two values, which were what got used in
    % different versions of the MScan software.
    if trimBorders
      margIm = mean(regImage{ch}, 1);
      LBWidth = find(margIm > 0, 1) - 1;
      RBWidth = imWidth - find(margIm > 0, 1, 'last');
      bWidthEst = min([LBWidth RBWidth]);
      [minBorderError, probableBorderI] = min(abs(validBorderWidths - bWidthEst));
      bWidth = validBorderWidths(probableBorderI);
      if minBorderError > 0
        warning('Possible border detection error: guess was %d, correcting to %d', bWidthEst, bWidth);
      end
      validPixels = [false(1, bWidth) true(1, imWidth - 2 * bWidth) false(1, bWidth)];
    end
    regImage{ch} = regImage{ch}(:, validPixels);

end
