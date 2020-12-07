function video_mergeVideos_ffmpeg(inputDirectory,inputExtension)
%% Merge an entire directory of behavior videos
% PAY ATTENTION TO THE CASE (UPPER OR LOWER) OF THE EXTENSION (THIS SOUNDS SILLY BUT IT MATTERS)

% fileDir= 'Z:\BehaviorVideo\mSM85\SpatialDisc\Session Data\mSM85_SpatialDisc_Jan28_2020_Session1';
iExt = inputExtension;
fileDir = inputDirectory;

fileListFace = dir(fullfile(fileDir,['*_1.' iExt ])); % specify the file extension to exclude directories
fileListMouth = dir(fullfile(fileDir,['*_2.' iExt])); % specify the file extension to exclude directories

if ~isempty(fileListFace)
    sepFileName = regexp(fileListFace(1).name,'_','split');
else
    sepFileName = regexp(fileListMouth(1).name,'_','split');
end
fileName = [sepFileName{1} '_' sepFileName{3} '_' sepFileName{4} '_' sepFileName{5}];

if sum(isspace(fileDir)) > 0
ffFileDir = regexp(fileDir,' ','split');
ffFileDir = [ffFileDir{1} '%20' ffFileDir{2}];
else 
    ffFileDir = fileDir;
end

for i = 1:length(fileListFace)
   fileListFace(i).faceText = ['file ' '''' fileDir filesep fileListFace(i).name ''''];
end
if ~isempty(fileListFace)
[fid,msg] = fopen(fullfile(fileDir,'faceText.txt'),'w');
assert(fid>=3,msg)
fprintf(fid,'%s\n',fileListFace.faceText)
fclose(fid);
end

fileListMouth = dir(fullfile(fileDir,['*_2.' iExt])); % specify the file extension to exclude directories
for i = 1:length(fileListMouth)
   fileListMouth(i).mouthText = ['file ' '''' fileDir filesep fileListMouth(i).name ''''];
end
if ~isempty(fileListMouth)
[fid,msg] = fopen(fullfile(fileDir,'mouthText.txt'),'w');
assert(fid>=3,msg);
fprintf(fid,'%s\n',fileListMouth.mouthText);
fclose(fid);
end

if iExt == 'mj2'
    commandFace = ['ffmpeg -f concat -safe 0 -i "' fileDir filesep 'faceText.txt" -c copy "' fileDir filesep fileName '_face.mj2"'];
    commandMouth = ['ffmpeg -f concat -safe 0 -i "' fileDir filesep 'mouthText.txt" -c copy "' fileDir filesep fileName '_mouth.mj2"'];
elseif iExt == 'mp4'
    commandFace = ['ffmpeg -f concat -safe 0 -i "' fileDir filesep 'faceText.txt" -c copy "' fileDir filesep fileName '_face.mp4"'];
    commandMouth = ['ffmpeg -f concat -safe 0 -i "' fileDir filesep 'mouthText.txt" -c copy "' fileDir filesep fileName '_mouth.mp4"'];
end
tic
[status,cmdout] = system(commandFace);
[status,cmdout] = system(commandMouth);
toc
