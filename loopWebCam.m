totalNumFrames = 300;
timeStamp = zeros(totalNumFrames,length(clock));
timeStamp1 = zeros(totalNumFrames,length(clock));

vidWriter = VideoWriter('frames.avi');
open(vidWriter);

for idx = 1:totalNumFrames
   % Acquire a single image.
   timeStamp(idx,:) = clock;
   rgbImage = snapshot(cam);
   timeStamp1(idx,:) = clock;
   writeVideo(vidWriter,rgbImage);
end

close(vidWriter);