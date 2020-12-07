clear all;

% matlab.video.read.UseHardwareAcceleration('off');
vidDir='Y:\xisun\videotest';
MJ2Ext='.mj2'; MP4Ext='.mp4';
fListMJ2=dir([vidDir filesep '*1' MJ2Ext]);
fListMP4=dir([vidDir filesep '*1' MP4Ext]);
convMetrics=struct;

 i=1;

if length(fListMJ2)==length(fListMP4)
    for i=1:length(fListMP4)
        disp(['Opening file #' num2str(i)]);
        rawData = squeeze(importdata([vidDir filesep fListMP4(i).name]));

        vMP4=VideoReader([vidDir filesep fListMP4(i).name]);
        fMP4=zeros(vMP4.Height,vMP4.Width,1,vMP4.NumFrames);
        fIdx=1;
        while(hasFrame(vMP4))
            fMP4(:,:,1,fIdx) = rgb2gray(readFrame(vMP4));
            fIdx=fIdx+1;
        end
        
        vMJ2=VideoReader([vidDir filesep fListMJ2(i).name]);
        fMJ2=zeros(vMJ2.Height,vMJ2.Width,1,vMJ2.NumFrames);
        fIdx=1;
        while(hasFrame(vMJ2))
            fMJ2(:,:,1,fIdx) = readFrame(vMJ2);
            fIdx=fIdx+1;
        end
        
        if i==1
            vidMP4=fMP4;
            vidMJ2=fMJ2;
        else
            vidMP4=cat(4,vidMP4,fMP4);
            vidMJ2=cat(4,vidMJ2,fMJ2);
        end
        
        convMetrics(i).pixelDifference = sum(abs(fMJ2(:)-fMP4(:)));
        convMetrics(i).fMJ2=fMJ2; convMetrics(i).fMP4=fMP4;
    end
    
else
    disp('The number of MJ2 vs MP4 files differs, review batch conver again')
end