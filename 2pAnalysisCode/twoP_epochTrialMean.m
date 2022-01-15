function eData = twoP_epochTrialMean(data,idxEpochs)
% Inputs:
% (1) data is a 3D matrix of cell (ROI) x frame x trial
% (2) idxEpochs epochs defines the frames start (column 1) and end (column 2) of each
% epoch (rows)

eData= zeros(size(data,1),size(idxEpochs,1),size(data,3));
for i = 1:size(idxEpochs,1)
    for j = 1:size(data,3)
        eData(:,i,j)= mean(data(:,idxEpochs(i,1):idxEpochs(i,2),j),2,'omitnan');
    end
end