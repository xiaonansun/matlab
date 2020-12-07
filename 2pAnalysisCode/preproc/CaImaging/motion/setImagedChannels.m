function channelsSaved = setImagedChannels(imageInfo)
% channelsSaved = setImagedChannels(imageInfo)
% Determins what scanning channels were saved.
%
% INPUT:
% imageInfo -- structure array of size 1x1, set by imfinfo(tifName); it
%        must be the imfinfo of only 1 frame (not the entire movie).


%%
if ~isfield('imageInfo', 'ImageDescription')
    warning('setImagedChannels:noImgDesc', 'imfinfo does not include ImageDescription; setting channelsSaved to 2!')
    warning('off','setImagedChannels:noImgDesc')
    channelsSaved = [1, 2];
    
else
    channelsSaved = [];
    if ~isempty(strfind(imageInfo.ImageDescription, 'Channel 1: Saved'))
        channelsSaved = [channelsSaved, 1];
    end
    
    if ~isempty(strfind(imageInfo.ImageDescription, 'Channel 2: Saved'))
        channelsSaved = [channelsSaved, 2];
    end
    
    if ~isempty(strfind(imageInfo.ImageDescription, 'Channel 3: Saved'))
        channelsSaved = [channelsSaved, 3];
    end
    
    if ~isempty(strfind(imageInfo.ImageDescription, 'Channel 4: Saved'))
        channelsSaved = [channelsSaved, 4];
    end
    
end


