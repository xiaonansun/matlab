function video = WebcamRecord(n, C)

save_dir='C:\Users\Anne\Documents\WebcamVideos';
clear cam;

cam=C;
cam.Resolution='640x480'; cam.FocusMode='manual'; cam.Focus=10; cam.Brightness=200; cam.Saturation=158; cam.Exposure=-5; cam.Contrast=50;

for i = 1:n
    rgbimage=snapshot(cam);
    video(:,:,i)=rgb2gray(rgbimage);
end
