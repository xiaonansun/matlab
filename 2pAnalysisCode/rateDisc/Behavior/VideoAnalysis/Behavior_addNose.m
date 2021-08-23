function SessionData = Behavior_addNose(fPath, SessionData, eyeCam)
% short code to assign snout and nose position to behavioral data.
if ~strcmpi(fPath(end),filesep)
    fPath = [fPath filesep];
end

if isempty(SessionData)
    noseSize = 60; %half-size of nose rectangle if reassigned
    SessionData.snoutPos(3) = 100; %half-size of nose rectangle if reassigned
    SessionData.eyePos(3) = 75; %half-size of nose rectangle if reassigned
else
    noseSize = 20; %half-size of nose rectangle if reassigned
end
    

if ~exist('eyeCam', 'var')
    eyeCam = 1; %camera that is viewing the face. Default is cam1.
end

%% check files and load raw data
movieFiles = dir([fPath '*Video*_' num2str(eyeCam) '.mp4']);

%% check if nose position is assigned or nose/snout position should be re-assigned
SessionData.eyePos(4) = eyeCam;
SessionData.snoutPos(4) = eyeCam;
SessionData.nosePos(3) = noseSize;
SessionData.nosePos(4) = eyeCam;

h1 = figure;
v = VideoReader([fPath movieFiles(1).name]);
temp = readFrame(v); imgSize = size(temp);clear v
imagesc(temp); axis image; colormap gray
text(round(imgSize(1) * 0.25), imgSize(2) - round(imgSize(2)/1.1),'click EYE','FontSize',30,'color','k')
[xPos,yPos] = ginput(1);
SessionData.eyePos(1) = round(xPos);
SessionData.eyePos(2) = round(yPos);

h2 = figure;
ax = subplot(1,3,1);
imagesc(temp);axis image; colormap gray
Frame = [SessionData.eyePos(1)-SessionData.eyePos(3) SessionData.eyePos(2)-SessionData.eyePos(3) SessionData.eyePos(3)*2 SessionData.eyePos(3)*2];
rectangle('Position',Frame,'linewidth',2,'edgecolor','k','parent',ax);
title('Eye position')

figure(h1);
v = VideoReader([fPath movieFiles(1).name]);
temp = readFrame(v);imgSize = size(temp);clear v
imagesc(temp);axis image; colormap gray
text(round(imgSize(1) * 0.25), imgSize(2) - round(imgSize(2)/1.1),'click SNOUT','FontSize',30,'color','b')
[xPos,yPos] = ginput(1);
SessionData.snoutPos(1) = round(xPos);
SessionData.snoutPos(2) = round(yPos);

figure(h2);
ax = subplot(1,3,2);
imagesc(temp);axis image; colormap gray
Frame = [SessionData.snoutPos(1)-SessionData.snoutPos(3) SessionData.snoutPos(2)-SessionData.snoutPos(3) SessionData.snoutPos(3)*2 SessionData.snoutPos(3)*2];
rectangle('Position',Frame,'linewidth',2,'edgecolor','b','parent',ax);
title('Snout position')

figure(h1);
v = VideoReader([fPath movieFiles(1).name]);
temp = readFrame(v);imgSize = size(temp); clear v
imagesc(temp);axis image; colormap gray
text(round(imgSize(1) * 0.25), imgSize(2) - round(imgSize(2)/1.1),'click NOSE','FontSize',30,'color','r')
[xPos,yPos] = ginput(1);
SessionData.nosePos(1) = round(xPos);
SessionData.nosePos(2) = round(yPos);

figure(h2);
ax = subplot(1,3,3);
imagesc(temp);axis image; colormap gray
Frame = [SessionData.nosePos(1)-SessionData.nosePos(3) SessionData.nosePos(2)-SessionData.nosePos(3) SessionData.nosePos(3)*2 SessionData.nosePos(3)*2];
rectangle('Position',Frame,'linewidth',2,'edgecolor','r','parent',ax);
title('Nose position');

close(h1); drawnow;