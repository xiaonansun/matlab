% This script outputs the framerates of the camera(s) attached to the computer.
% Cameras are indexed form 0, therefore the first element of the array refers to the camera with camera ID of 0

commandStr = 'python "C:\Users\Anne\Dropbox\Users\Richard\python_scripts\GetCameraFrameRate.py"';
[status, commandOut] = system(commandStr);
camera_framerates = round(str2num(commandOut));