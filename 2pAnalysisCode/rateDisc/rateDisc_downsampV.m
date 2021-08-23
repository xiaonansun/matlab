function newV = rateDisc_downsampV(V,b)
% downsample temporal component V of the form (dims x frames x trials) to a
% new V. b determines the downsampling factor (1/b).

newV = NaN(size(V,1), floor(size(V,2)*(1/b)), size(V,3), 'single');
for iTrials = 1:size(V,3)
    
    cData = double(squeeze(V(:,:,iTrials)))';
    cIdx = ~isnan(cData(:,1));
    dShift = find(cIdx,1); %shift data if trial is starting with NaNs
    
    if dShift > 1
        cIdx(1 : dShift + rem(dShift,b)) = false;
        dShift = (find(cIdx,1)-1) / b + 1;
    end
    
    newData = resample(cData(cIdx,:) - mean(cData(cIdx,:)),1,b) + mean(cData(cIdx,:)); %downsample
    newV(:,dShift - 1 + (1:size(newData,1)), iTrials) = newData';
    
%     plot(0:size(newV,2)-1,newV(1,:,iTrials)); hold on;
%     plot(0:0.5:size(V,2)/2-0.5,V(1,:,iTrials)); hold off; pause;
end
end