animal = 'Plex01';
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
cropRange = 150; %center range of the image used for for second alignment
modelPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\models\190217_213418-n=103\final_model.h5'; %leap model for initial alignment

%% get recordings
fPath = [cPath animal filesep 'SpatialDisc' filesep];
recs = ls(fPath);
recs = recs(~ismember(recs(:,1), '.'), :);

load([fPath filesep strtrim(recs(1,:)) filesep 'firstAlign.mat'], 'firstRotate', 'firstTranslate') %load initial alignment for reference
load([fPath strtrim(recs(1,:)) filesep 'blueAvg.mat']);
params.angleD = firstRotate;
params.scaleConst = 1;
params.tC = firstTranslate;
alignPic = alignAllenTransIm(blueAvg, params);
alignRec = strtrim(recs(1,:));

%% for new recordings run pre-alignment first
if strcmpi(type, 'new')
    
    recs = ls(fPath);
    recs = recs(~ismember(recs(:,1), '.'), :);
    snapMovie = NaN(size(blueAvg,1), size(blueAvg,2), size(recs,1), 'single');
    cIdx = false(1, size(recs,1));
    
    for iRecs = 1 : size(recs,1)
        try
            if exist([fPath strtrim(recs(iRecs,:)) filesep 'opts2.mat'], 'file') == 2 %check for alignment
                cIdx(iRecs) = true;
            else
                load([fPath strtrim(recs(iRecs,:)) filesep 'blueAvg.mat']);
                snapMovie(:, :, iRecs) = blueAvg;
            end
        catch
            cIdx(iRecs) = true;
        end
    end
    recs(cIdx, :) = [];
    snapMovie(:, :, cIdx) = [];
    clear cIdx
    
    if isempty(snapMovie) %no new recordings
        error('No new recordings found. Use type = "exist" to re-align existing sessions.')
    end
    
    %% align images using leap network
    tempMovie = arrayResize(snapMovie, 2);
    tempMovie(272, :, :, :) = 0;
    tempMovie = reshape(tempMovie, size(tempMovie,1), size(tempMovie,2), 1, []);
    
    preds = predict_box(uint8((tempMovie(:,:,1,:)./ 65535) .* 255), modelPath);
    lamb2breg = squeeze(preds.positions_pred(4, :, :) - preds.positions_pred(5, :, :));
    lamb2front = squeeze(preds.positions_pred(2, :, :) - preds.positions_pred(5, :, :));
    bregfront = squeeze(preds.positions_pred(2, :, :) - preds.positions_pred(4, :, :));
    
    cRotate = atan2(lamb2breg(2,:), lamb2breg(1,:)) .* 180 ./ pi .* -1;
    cRotate(2,:) = atan2(lamb2front(2,:), lamb2front(1,:)) .* 180 ./ pi .* -1;
    cRotate(3,:) = atan2(bregfront(2,:), bregfront(1,:)) .* 180 ./ pi .* -1;
    cRotate = reshape(cRotate, size(cRotate,1),[]);
    cRotate = 90 - mean(cRotate); %mean angle of center line
    
    bregmaPos = reshape(preds.positions_pred(4, :, :), 2, size(snapMovie,3)); %mean position of bregma
    bregmaPos = squeeze(bregmaPos) * 2;
    
    lambdaPos = reshape(preds.positions_pred(5, :, :), 2, size(snapMovie,3)); %mean position of lambda
    lambdaPos = squeeze(lambdaPos) * 2;
    
    frontPos = reshape(preds.positions_pred(2, :, :), 2, size(snapMovie,3)); %mean position of front center
    frontPos = squeeze(frontPos) * 2;
    
    %% make first aligned movie
    alignMovie = NaN(size(snapMovie,1),size(snapMovie,2),size(snapMovie,3), 'single');
    cTranslate = NaN(2, size(snapMovie,3));
    for iFrames = 1 : length(cRotate)
        
        temp = snapMovie(:,:,iFrames);
        dSize = size(temp);
        [xx,yy] = meshgrid(1:dSize(2),1:dSize(1)); %isolate index for selected area
        mask = false(dSize(1:2));
        mask = mask | hypot(xx - bregmaPos(1, iFrames), yy - bregmaPos(2, iFrames)) <= 2;
        temp(mask) = Inf;
        
        mask = false(dSize(1:2));
        mask = mask | hypot(xx - lambdaPos(1, iFrames), yy - lambdaPos(2, iFrames)) <= 2;
        mask = mask | hypot(xx - frontPos(1, iFrames), yy - frontPos(2, iFrames)) <= 2;
        temp(mask) = NaN;
        
        temp = imrotate(temp, cRotate(iFrames), 'crop');
        xCenter = round(mean(ceil(find(isnan(temp) | isinf(temp))/size(temp,1) ))); %find center line
        yCenter = rem(find(isinf(temp),1), size(temp,1)); %find bregma
        
        cTranslate(:, iFrames) = [round(size(temp,2)/2) - xCenter, round(size(temp,1)/2) - yCenter];
        alignMovie(:,:,iFrames) = imrotate(snapMovie(:,:,iFrames), cRotate(iFrames), 'crop');
        alignMovie(:,:,iFrames) = imtranslate(alignMovie(:,:,iFrames), cTranslate(:, iFrames));
    end
    
    %% align to first frame one more time to reduce variability
    tempMovie = gpuArray(alignMovie);
    alignPic = gpuArray(alignPic);
    
    options.TolX = 0.0001;
    for iFrames = 1 : size(tempMovie,3)
        [nRotate(iFrames), nError(iFrames)] = fminsearch(@(u) Widefield_alignImage(alignPic, tempMovie(:,:,iFrames), u, cropRange), 1, options);
        [~, nTranslate(:,iFrames), tempMovie(:, :, iFrames)] = Widefield_alignImage(alignPic, tempMovie(:, :, iFrames), nRotate(iFrames), 0);
    end
    finalRotate = cRotate + nRotate;
    finalTranslate = round(cTranslate + flipud(nTranslate));
    
    alignPic = gather(alignPic);
    clear tempMovie
    
else
    % get snapshot based on last run of alignment code
    snapMovie = NaN(size(blueAvg,1), size(blueAvg,2), size(recs,1), 'single');
    useIdx = true(1, size(recs,1));
    for iRecs = 1 : size(recs,1)
        try
            load([fPath strtrim(recs(iRecs,:)) filesep 'blueAvg.mat']);
            snapMovie(:, :, iRecs) = blueAvg;
        catch
            % dont use this recording - throw warning
            disp(['!!! Couldnt use ' strtrim(recs(iRecs,:)) ' - removed from alignment.!!!']);
            useIdx(iRecs) = false;
        end
    end
    recs = recs(useIdx,:);
    snapMovie = snapMovie(:,:,useIdx);
    finalRotate = finalRotate(useIdx); 
    finalTranslate = finalTranslate(:,useIdx); 
    finalScale = finalScale(useIdx);
end

%% align images and run inspection tool
firstRotate = finalRotate;
firstTranslate = finalTranslate;
firstScale = finalScale;
for iFrames = 1 :size(snapMovie,3)
    params.angleD = firstRotate(iFrames);
    params.scaleConst = firstScale(iFrames);
    params.tC = firstTranslate(:,iFrames);
    temp(:,:,iFrames) = alignAllenTransIm(snapMovie(:,:,iFrames), params);
end
temp(isnan(temp)) = 0; %don't keep NaNs
BrainAligner(alignPic, temp); % use this to create appRotate and appShift

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run this code when things look fine
finalTranslate = firstTranslate + appShift;
finalRotate = firstRotate + appRotate;
finalScale = firstScale + appScale - 1;

for iFrames = 1 : size(snapMovie,3)
    params.angleD = finalRotate(iFrames);
    params.scaleConst = finalScale(iFrames);
    params.tC = finalTranslate(:, iFrames);
    temp2(:,:,iFrames) = alignAllenTransIm(snapMovie(:,:,iFrames), params);
end
compareMovie(temp2); %double-check alignment

%% save into individual folders
for iRecs = 1 : size(recs,1)
    firstRotate = finalRotate(iRecs);
    firstTranslate = finalTranslate(:, iRecs);    
    save([fPath strtrim(recs(iRecs,:)) filesep 'firstAlign.mat'], 'firstTranslate', 'firstRotate')
end

%% load allen alignment in reference recording
load([fPath alignRec filesep 'opts2.mat'])
load([fPath alignRec filesep 'blueAvg.mat'])
load([fPath alignRec filesep 'firstAlign.mat'])

tC = opts.transParams.tC - firstTranslate;
angleD = opts.transParams.angleD - firstRotate;
allenMovie(:,:,1) = alignAllenTransIm(blueAvg, opts.transParams);

Cnt = 0;
for iRecs = find(~contains(cellstr(recs),alignRec))'
    load([fPath alignRec filesep 'opts2.mat'])
    load([fPath strtrim(recs(iRecs,:)) filesep 'firstAlign.mat'])
    load([fPath strtrim(recs(iRecs,:)) filesep 'blueAvg.mat'])
    
    opts.transParams.tC = tC + firstTranslate;
    opts.transParams.angleD = angleD + firstRotate;
    
    save([fPath strtrim(recs(iRecs,:)) filesep 'opts2.mat'], 'opts')
    Cnt = Cnt + 1;
    allenMovie(:,:,Cnt) = alignAllenTransIm(blueAvg, opts.transParams);
end
compareMovie(allenMovie); %double check allen alignment

%% compare with previous overview alignment file and update


% SOMETHING IS WRONG HERE !!
old = load([cPath animal filesep 'firstAlignOpts.mat'], 'finalRotate', 'finalTranslate', 'recs');
for iRecs = 1 : size(recs,1)
    tIdx = contains(cellstr(old.recs), recs(iRecs,:));
    if any(tIdx)
        if ~strcmpi(old.recs(tIdx,:), recs(iRecs,:))
            error('Just to be sure :)');
        end
        old.finalRotate(tIdx) = finalRotate(iRecs);
        old.finalTranslate(:, tIdx) = finalTranslate(:, iRecs);
    else
        old.recs(end+1,:) = recs(iRecs,:); %make new entry
        old.finalRotate(end+1) = finalRotate(iRecs); %make new entry
        old.finalTranslate(:, end+1) = finalTranslate(:, iRecs); %make new entry
    end
end
recs = old.recs;
finalRotate = old.finalRotate;
finalTranslate = old.finalTranslate;
save([cPath animal filesep 'firstAlignOpts.mat'], 'finalRotate', 'finalTranslate', 'recs') %save into base folder
