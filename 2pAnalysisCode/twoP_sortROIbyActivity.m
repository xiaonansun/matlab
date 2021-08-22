function [oMatrix1,oMatrix2,oMatrix12,oMatrix21]=twoP_sortROIbyActivity(iMatrix1,iMatrix2,filterType,filterWindow)
% This function sorts the PSTH for multiple ROIs based on the timing of
% peak activity
% The input is a 2-D matrix and the output is a 2-D matrix of identical
% dimension.
% iMatrix: the input 2-D matrix where rows are ROIs and columns represent
% time
% filterType: type of smoothing filter to apply. The default is
% the Savitzky-Golay (polynomial) filter.
% filterWindow: the smoothing filter window. Default is 10 data points.

if ~exist('filterType','var')
    filterType = 'sgolay';
end

if ~exist('filterWindow','var') 
    filterWindow = 10;
end

if ~exist('iMatrix2','var') || isempty(iMatrix2)
    [~,maxValIdx1]= max(smoothdata(iMatrix1,2,filterType,filterWindow),[],2); maxValIdx1= sortrows([maxValIdx1 (1:1:length(maxValIdx1))']);
    oMatrix1 = iMatrix1(maxValIdx1(:,2),:);
    oMatrix2 = zeros(size(oMatrix1,1),size(oMatrix1,2));
else
    [~,maxValIdx1]= max(smoothdata(iMatrix1,2,filterType,filterWindow),[],2); maxValIdx1= sortrows([maxValIdx1 (1:1:length(maxValIdx1))']);
    oMatrix1 = iMatrix1(maxValIdx1(:,2),:);
    [~,maxValIdx2]= max(smoothdata(iMatrix2,2,filterType,filterWindow),[],2); maxValIdx2= sortrows([maxValIdx2 (1:1:length(maxValIdx2))']);
    oMatrix2 = iMatrix2(maxValIdx2(:,2),:);
    oMatrix12 = iMatrix1(maxValIdx2(:,2),:);
    oMatrix21 = iMatrix2(maxValIdx1(:,2),:);
end

end