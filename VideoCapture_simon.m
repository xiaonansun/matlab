adapter = 'winvideo';
camSet = 'YUY2_640x480';
vid=videoinput(adapter,1,camSet);
vid_info=imaqhwinfo(vid);

% BpodSystem.PluginObjects.Webcams = videoinput(adapter,1,camSet);
