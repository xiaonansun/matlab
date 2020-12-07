clear all; clc;

%% OPEN PREVIEW
deviceIdx = 2;
formatIdx = 4; % Cell index of the format

adaptorInfo = imaqhwinfo;
deviceID = imaqhwinfo(adaptorInfo.InstalledAdaptors{deviceIdx});
deviceInfo = imaqhwinfo(adaptorInfo.InstalledAdaptors{deviceIdx},1);
vid = videoinput(adaptorInfo.InstalledAdaptors{deviceIdx},1,deviceInfo.SupportedFormats{formatIdx});

%% ACQUIRE IMAGE DATA
% Create video input object. 
vid = videoinput(adaptorInfo.InstalledAdaptors{deviceIdx},1,deviceInfo.SupportedFormats{formatIdx})

% Set video input object properties for this application.
vid.FramesPerTrigger = 2000;
vid.LoggingMode = 'disk'; 

logfile = VideoWriter('logfile.mp4', 'Grayscale AVI');
vid.DiskLogger = logfile;

start(vid);
wait(vid);

%% PREVIEW VIDEO

preview(vid);

%% CLOSE PREVIEW

closepreview(vid); 