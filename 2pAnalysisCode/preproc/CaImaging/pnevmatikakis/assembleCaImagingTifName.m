function fName = assembleCaImagingTifName(nums, oldTifName, badFrames)
% fName = assembleCaImagingTifName(nums [, oldTifName] [, badFrames])
% 
% Take the nums output from parseCaImagingTifName() and turn it back into a
% filename.
%
% If oldTifName is specified and true, use the old-style tif naming
% convention (2-digit minor number). Default false.
%
% If badFrames is specified and true, the filename will end with
% _badFrames.mat . This overrides the _MCM option below. Note that
% badFrames is no longer used in the current extraction pipeline.
%
% If nums(4) is non-NaN, the filename will end with _ch#_MCM.TIF . If it is
% NaN, it will not include the channel number or MCM.


%% Optional oldTimeName argument

if ~exist('oldTifName', 'var') || isempty(oldTifName)
    oldTifName = false;
end


%% Tif name stem

if ~isnan(nums(4))
    if oldTifName==1
        fName = sprintf('%06d_%03d_%02d_ch%d', nums(1), nums(2), nums(3), nums(4));
    elseif oldTifName==2
        fName = sprintf('%06d_%03d_%01d_ch%d', nums(1), nums(2), nums(3), nums(4));
    else
        fName = sprintf('%06d_%03d_%03d_ch%d', nums(1), nums(2), nums(3), nums(4));
    end
else
    if oldTifName==1
        fName = sprintf('%06d_%03d_%02d', nums(1), nums(2), nums(3));
    elseif oldTifName==2
        fName = sprintf('%06d_%03d_%01d', nums(1), nums(2), nums(3));
    else
        fName = sprintf('%06d_%03d_%03d', nums(1), nums(2), nums(3));
    end
end


%% End of name (badFrames or MCM etc.)

if exist('badFrames', 'var') && badFrames
  fName = [fName '_badFrames.mat'];
elseif ~isnan(nums(4)) % nums(5) == 1
  fName = [fName '_MCM.TIF'];
else
  fName = [fName '.TIF'];
end
