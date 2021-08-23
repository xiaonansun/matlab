function Frame = Behavior_rebuildSegmentFrame(fPath, eyeCam)

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end
load([fPath '\segInd'  num2str(eyeCam) '.mat'], 'ind')

Frame = NaN(240,320);
for iSegs = 1 : size(ind,2)
    
    load([fPath 'SVD_Cam' int2str(eyeCam) '-Seg' int2str(iSegs) '.mat'],'U','V'); %load current segment
    Frame(ind(:,iSegs)) = V(1,:) * U;
    
end
Frame = imresize(Frame,2);