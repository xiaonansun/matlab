function res = lightsheet_averageStack(sampleID, nImg, zProjectMethod)
% This function averages a specified number of consecutive images (nImg)
% and forms a smaller tif stack

fileExt = '.tif';
j2kExt = '.j2k';

opts.big = true;
opts.append = true;
opts.message = false;

baseDir = 'F:\LightsheetBatch2020-10-18';
imgDir = fullfile(baseDir,sampleID); 
imgDirContent = dir(fullfile(imgDir,['*' j2kExt])); imgDirContent = natsortfiles({imgDirContent.name});
zprjDir = fullfile(imgDir,'zProjectedStack');

if ~exist(zprjDir,'dir')
    mkdir(zprjDir);
else
    zprjDirCont = dir(fullfile(zprjDir,'*substack.tif'));
    if ~isempty(zprjDirCont)
       for i = 1:length(zprjDirCont) 
           zprjFilePath = fullfile(zprjDirCont(i).folder,zprjDirCont(i).name);
           delete(zprjFilePath);
       end
    end
end

if ~exist('zProjectMethod','var') || isempty(zProjectMethod)
    zProjectMethod = 'maximum';
end
disp(['The ' zProjectMethod ' from every ' num2str(nImg) ' consecutive images will be computed to generate a new stack with ' num2str(ceil(length(imgDirContent)/nImg)) ' slices.']);

tempImgStack = imread(fullfile(imgDir,imgDirContent{1}));
nRows = size(tempImgStack,1); nCols = size(tempImgStack,2); imgClass = whos('tempImgStack'); imgClass = imgClass.class;
clear tempImgStack;

if rem(length(imgDirContent),nImg)
    idxA = 1:nImg:length(imgDirContent); idxA = idxA(1:end);
    idxB = 0:nImg:length(imgDirContent); idxB = [idxB(2:end) length(imgDirContent)];
else
    idxA = 1:nImg:length(imgDirContent);
    idxB = 0:nImg:length(imgDirContent); idxB = idxB(2:end);
end

if length(idxA) ~= length(idxB)
    disp('Error: there is mismatch in the number of image indices! Check idxA and idx B. Script terminating')
    return
end

% C = repmat(1:1:nImg,length(idxA),1);
D = reshape(1:1:ceil(length(imgDirContent)/nImg)*nImg,[nImg,ceil(length(imgDirContent)/nImg)])';
% zprjImg = zeros(nRows,nCols,ceil(length(imgDirContent)/nImg),imgClass);
for i = 1:length(idxA)
    tempImgStack = zeros(nRows,nCols,length(idxA(i):1:idxB(i)),imgClass);
    disp(['Loading images ' num2str(idxA(i)) ' to ' num2str(idxB(i))]);
%     sC = C(i,:);
    sD = D(i,:);
    parfor j = 1:length(idxA(i):1:idxB(i))
        k = sD(j);
        fprintf(['Loading image ' imgDirContent{k} '...']);
        tempImgStack(:,:,j) = imread(fullfile(imgDir,imgDirContent{k}));
        fprintf('Done! \n');
    end
    disp(['Averaging image ' num2str(i)]);
    zprjImg = max(tempImgStack,[],3);
    disp(['Saving z-projected image #' num2str(i)]);
    res = saveastiff(zprjImg,fullfile(zprjDir, [sampleID '_zProject_' num2str(nImg) '_' zProjectMethod '_substack.tif']),opts);
    clear zprjImg tempImgStack;
end

% imwrite(zprjImg, fullfile(dsImgDir, [sampleID '_zProject_' num2str(nImg) '_' zProjectMethod '_substack.tif']));
