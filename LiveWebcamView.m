%%Example WEBCAM custom preview window for MATLAB R2017a
clear all;

%%List connected webcams
webcamlist

%%Connect to webcam
c = webcam(1);

%%Setup preview window
fig = figure('NumberTitle', 'off', 'MenuBar', 'none');
fig.Name = 'My Camera';
ax = axes(fig);
frame = snapshot(c);
im = image(ax, zeros(size(frame), 'uint8'));
axis(ax, 'image');

%%Start preview
preview(c, im)
setappdata(fig, 'cam', c);
fig.CloseRequestFcn = @closePreviewWindow_Callback;

%%Recording
save_dir='C:\Users\Anne\Documents\MATLAB';
video_filename='webcam_video.dat';
timestamp_filename = 'webcam_timestamp.dat';
video_dataID=fopen([save_dir '\' video_filename],'w');
timestampID=fopen([save_dir '\' timestamp_filename],'w');
n=1;
while n<=90
    [rgb_image timestamp]=snapshot(c);
    gray_image=uint8(rgb2gray(rgb_image));
    fwrite(video_dataID,gray_image,'uint8');
    fwrite(timestampID,timestamp,'double');
    n=n+1;
end

%%Local functions
function closePreviewWindow_Callback(obj, ~)
c = getappdata(obj, 'cam');
closePreview(c)
delete(obj)
end