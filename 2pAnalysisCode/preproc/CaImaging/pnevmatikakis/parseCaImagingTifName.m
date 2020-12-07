function [nums, valid, oldTifName, renameTifAddMinor] = parseCaImagingTifName(fName)
% [nums, valid] = parseCaImagingTifName(fName)
% 
% This function parses a .TIF filename fName of the form YYMMDD_mmm_nn.TIF
% or YYMMDD_mmm_nn_ch#_MCM.TIF, where YYMMDD is a 6-digit date, mmm is the
% three-digit "major" number, and nn is the 2-digit "minor" number. ch#
% indicates 1-digit channel number which is only present in
% motion-corrected TIF file names. If _MCM is present, this will be noted
% in the output.

% nums   -- numeric array containing [date, major, minor, channel, hadMCM]
%           where hadMCM is a boolean indicating whether the filename ended
%           in _MCM.TIF or just .TIF. channel will be NaN for tif file
%           names not including a channel number (ie the case with raw tif
%           files before motion correction).
% valid  -- boolean indicating whether the provided name had the correct
%           form. If not, nums will be a 1 x 5 of NaNs
% oldTifName -- boolean, indicating whether tif files are named in the old
%           way, ie tif minor is XX (instead of XXX). 
%
% To turn these numbers back into a filename, use
% assembleCaImagingTifName()


% Note: this could have been implemented more simply with sscanf, but that
% wouldn't allow as much error checking.


%% Set default return values, in case invalid name
nums = NaN(1, 4);
valid = false;
oldTifName = false;

%% Pull off path if present, filename stem, extension
[~, stem, ext] = fileparts(fName);


%% Tokenize filename using '_' as separator
tokens = simpleTokenize(stem, '_');


%% FN: in MView's automatic conversion, tif files won't have a tif minor sufix if total frame numbers of the mdf file < 4089.
% Here we assign tif minor 001 to them so the rest of the function goes
% through, and flag them so they can be renamed outside this function.

if length(tokens)==2 && ...
    (strcmpi(ext, '.tif') || strcmpi(ext, '.tiff') || strcmpi(ext, '.bz2')) && ...
    length(tokens{1}) == 6 && length(tokens{2}) == 3
    
    % Add 001 as tif minor for consistency in file name with all other tif files.
    tokens{3} = '001';
    fprintf('_001 added as tif minor to file %s\n', fName)
    renameTifAddMinor = 1; % flag for changing file name
else
    renameTifAddMinor = 0;
end


%% Initial check for validity

% Check: extension, number of tokens, lengths of tokens
if ~(strcmpi(ext, '.tif') || strcmpi(ext, '.tiff') || strcmpi(ext, '.bz2')) ...
    || ~ismember(length(tokens), [2 3 4 5]) ... 
    || length(tokens{1}) ~= 6 || length(tokens{2}) ~= 3 ...
    || ~ismember(length(tokens{3}), [1 2 3]) % FN changed this, since in the new version of mview tif files are named "151006_001_001.TIF"
  % Invalid name, return
  return;
end

% FN commented the following, since raw tif files won't have _ch in their
% names.
%{
if ~(strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')) || ...
    ~ismember(length(tokens), [4 5]) || ...
    length(tokens{1}) ~= 6 || length(tokens{2}) ~= 3 || ...
    length(tokens{3}) ~= 2 || length(tokens{4}) ~= 3
  % Invalid name, return
  return;
end
%}


%% Set oldTifName to 1, if tif file names are in the old format (ie minor is XX instead of XXX)
if length(tokens{3})==2
    oldTifName = 1;
elseif length(tokens{3})==1 % FN: one-digit tif minor exists for some old tif files (recorded around Aug 2015)
    oldTifName = 2;
else
    oldTifName = 0;
end

%% Convert date, major/minor numbers, channel to numeric, more checks

theDate = str2double(tokens{1});
major = str2double(tokens{2});
minor = str2double(tokens{3});
if length(tokens)>3
    channel = str2double(tokens{4}(3));
else
    channel = NaN;
end

if isnan(theDate) || isnan(major) || isnan(minor) % || isnan(channel) || ...
    % ~strcmp(tokens{4}(1:2), 'ch') 
  % Invalid name, return
  return;
end


%% Check for _MCM --> FN: we don't need this anymore, since MCM tif files have "_ch" in their name, hence already distinguishable from non MCM files.
%{
if length(tokens) == 5 && strcmpi(tokens{5}, 'MCM') 
  hadMCM = 1;
else
  hadMCM = 0;
end
%}

%% Assemble return values

valid = true;
nums = [theDate, major, minor, channel];
% nums = [theDate, major, minor, channel, hadMCM];


