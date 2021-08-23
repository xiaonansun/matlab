function h2 = Behavior_checkNose(pic, SessionData, eyeCam)
% short code to check snout and nose position as defined in SessionData.

if ~exist('eyeCam', 'var')
    eyeCam = 1; %camera that is viewing the face. Default is cam1.
end

%% check if nose position is assigned or nose/snout position should be re-assigned
h2 = figure;
ax = subplot(1,3,1);
imagesc(pic);axis image; colormap gray
Frame = [SessionData.eyePos(1)-SessionData.eyePos(3) SessionData.eyePos(2)-SessionData.eyePos(3) SessionData.eyePos(3)*2 SessionData.eyePos(3)*2];
rectangle('Position',Frame,'linewidth',2,'edgecolor','k','parent',ax);
title('Eye position')

figure(h2);
ax = subplot(1,3,2);
imagesc(pic);axis image; colormap gray
Frame = [SessionData.snoutPos(1)-SessionData.snoutPos(3) SessionData.snoutPos(2)-SessionData.snoutPos(3) SessionData.snoutPos(3)*2 SessionData.snoutPos(3)*2];
rectangle('Position',Frame,'linewidth',2,'edgecolor','b','parent',ax);
title('Snout position')

if isfield(SessionData, 'nosePos')
figure(h2);
ax = subplot(1,3,3);
imagesc(pic);axis image; colormap gray
Frame = [SessionData.nosePos(1)-SessionData.nosePos(3) SessionData.nosePos(2)-SessionData.nosePos(3) SessionData.nosePos(3)*2 SessionData.nosePos(3)*2];
rectangle('Position',Frame,'linewidth',2,'edgecolor','r','parent',ax);
title('Nose position');
end
drawnow;