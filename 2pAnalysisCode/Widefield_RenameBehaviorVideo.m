function Widefield_RenameBehaviorVideo(firstPath,secondPath,offset)
% code to rename behavioral video files if a recording was interrupted.
% This is to ensure continous trial enumeration where videos in secondPath
% use same file names as in firstPath and start with trialNr given by offset.

% firstPath = 'U:\space_managed_data\BehaviorVideo\mSM30\SpatialDisc\Session Data\mSM30_SpatialDisc_Aug10_2017_Session1';
% secondPath = 'U:\space_managed_data\BehaviorVideo\mSM30\SpatialDisc\Session Data\mSM30_SpatialDisc_Aug10_2017_Session2';

[~, b] = fileparts(firstPath);
temp = findstr(b,'_');
cAnimal = b(1:temp(1)-1);

orgVids = dir([firstPath filesep cAnimal '*mj2']);
orgVids = orgVids(1).name;
ind = strfind(orgVids,'_');
orgVids(ind(end-2):end) = [];

nextVids = dir([secondPath filesep cAnimal '*mj2']);
nextVids = nextVids(1).name;
ind = strfind(nextVids,'_');
nextVids(ind(end-2):end) = [];

vidFiles = dir([secondPath filesep cAnimal '*1.mj2']);
frameFiles = dir([secondPath filesep cAnimal '*1.mat']);

for iFiles = 1 : length(vidFiles)
    
    temp = strrep(vidFiles(iFiles).name, nextVids, orgVids);
    temp = textscan(vidFiles(iFiles).name,'%s%s%s%s%s%s%s%s','delimiter','_');
    cFile = [secondPath '\' orgVids '_Video_' num2str(iFiles + offset, '%04i') '_' temp{end}{1}]; %current video file
    movefile([secondPath '\' vidFiles(iFiles).name],cFile);
    
    temp = strrep(frameFiles(iFiles).name, nextVids, orgVids);
    temp = textscan(frameFiles(iFiles).name,'%s%s%s%s%s%s%s%s','delimiter','_');
    cFile = [secondPath '\' orgVids '_frameTimes_' num2str(iFiles + offset, '%04i') '_' temp{end}{1}]; %current video file
    movefile([secondPath '\' frameFiles(iFiles).name],cFile);
    
end

vidFiles = dir([secondPath filesep cAnimal '*2.mj2']);
frameFiles = dir([secondPath filesep cAnimal '*2.mat']);

for iFiles = 1 : length(vidFiles)
    
    temp = strrep(vidFiles(iFiles).name, nextVids, orgVids);
    temp = textscan(vidFiles(iFiles).name,'%s%s%s%s%s%s%s%s','delimiter','_');
    cFile = [secondPath '\' orgVids '_Video_' num2str(iFiles + offset, '%04i') '_' temp{end}{1}]; %current video file
    movefile([secondPath '\' vidFiles(iFiles).name],cFile);
    
    temp = strrep(frameFiles(iFiles).name, nextVids, orgVids);
    temp = textscan(frameFiles(iFiles).name,'%s%s%s%s%s%s%s%s','delimiter','_');
    cFile = [secondPath '\' orgVids '_frameTimes_' num2str(iFiles + offset, '%04i') '_' temp{end}{1}]; %current video file
    movefile([secondPath '\' frameFiles(iFiles).name],cFile);
    
end
