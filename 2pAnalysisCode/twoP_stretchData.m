function newVc = twoP_stretchData(Vc)
% needs to be fixed, the output is a matrix of zeros!!!
% This function stretches neural data, if needed, such that the duration of a given epoch matches across all trials 
S=twoP_settings;
segFrames = S.segFrames;

msPerFrame=32.3638;
frameRate = 1000/msPerFrame; % Frame rate of imaging
sRate = frameRate; 
% segFrames = cumsum(floor(segIdx * sRate)); % max nr of frames per segment
newVc = zeros(size(Vc,1),size(Vc,2),size(Vc,3));

for j = 1:size(Vc,3)
handle = Vc(:,1:segFrames(2),j); 
stim = Vc(:,segFrames(2)+1:segFrames(3),j);
delay = Vc(:,segFrames(3)+1:segFrames(4),j); 
response = Vc(:,segFrames(4)+1:segFrames(end),j); 

maxH = size(handle,2);
for i = 1:size(handle,1)
    d = sum(~isnan(handle(i,:)));
    x=1:d;
    v=handle(i,x);
    xq = linspace(1,d,maxH);
    intH(i,:)=interp1(x,v,xq);
end

maxS = size(stim,2);
for i = 1:size(stim,1)
    d = sum(~isnan(stim(i,:)));
    x=1:d;
    v=stim(i,x);
    xq = linspace(1,d,maxS);
    intS(i,:)=interp1(x,v,xq);
end

maxD = size(delay,2);
for i = 1:size(delay,1)
    d = sum(~isnan(delay(i,:)));
    if d < 2
        x=1:2;
        v=[delay(i,1) delay(i,1)];
    else
        x=1:d;
        v=delay(i,x);
    end
    xq = linspace(1,d,maxD);
    intD(i,:)=interp1(x,v,xq);
end
newVC(:,:,j) = [intH intS intD response];

end


%%
% histogram(sum(~isnan(handle),2)*data.msPerFrame);
% histogram(sum(~isnan(stim),2)*data.msPerFrame);
% histogram(sum(~isnan(delay),2)*data.msPerFrame,100);
% histogram(sum(~isnan(response),2)*data.msPerFrame,100);