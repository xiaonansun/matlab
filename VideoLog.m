vid = videoinput('winvideo', 1, 'YUY2_640x480');
src = getselectedsource(vid);
vid.FramesPerTrigger = Inf;
vid.ReturnedColorspace = 'grayscale';
src.BacklightCompensation = 'on';
src.ExposureMode = 'manual';
src.FocusMode = 'manual';
src.Exposure = -4;
src.Brightness = 162;
src.Contrast = 129;
src.Focus = 20;
src.Gain=255;
src.Saturation = 0;
src.WhiteBalance = 6000;
src.WhiteBalanceMode = 'manual';1 provided. No such port open. Maybe you closed it beforehand?
Error in function Close: 	Usage error
Invalid port handle -1 provided. 

%% previews vid object
preview(vid);

%%
vid.LoggingMode = 'disk';
diskLogger = VideoWriter('C:\Users\Anne\Documents\test_video_0001.avi', 'Uncompressed AVI');
vid.DiskLogger = diskLogger;
triggerconfig(vid, 'manual');

%% starts video
start(vid);

%% triggers acquisition
trigger(vid);

%% stops video acquisition
stop(vid);
