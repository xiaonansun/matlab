animal = 'CSP22';
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
cropRange = 150; %center range of the image used for for second alignment
modelPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\models\190217_213418-n=103\final_model.h5'; %leap model for initial alignment

%% get alignment pic - this should be created by aligning one recording to the allen CCF
load([cPath animal filesep 'alignPic.mat'])
load([cPath animal filesep 'alignOpts.mat'])
alignPic = alignAllenTransIm(blueAvg, opts.transParams);
alignOpts = opts;

%% get all recordings
fPath = [cPath animal filesep 'SpatialDisc' filesep];
recs = rateDisc_getRecsForAnimal(animal);
recs = cat(1,recs(:).name);

%% check for alignment data and create a preliminary opts file otherwise
rawMovie = NaN(size(alignPic,1), size(alignPic,2), size(recs,1), 'single'); %raw vessel images
snapMovie = NaN(size(alignPic,1), size(alignPic,2), size(recs,1), 'single'); %first snap movie
appShift = NaN(2, size(recs,1)); %translation values
appRotate = NaN(1, size(recs,1)); %rotation values
appScale = NaN(1, size(recs,1)); %scaling values

rejIdx = false(1, size(recs,1)); %index for broken recordings
for iRecs = 1 : size(recs,1)
    try
        cPath = [fPath strtrim(recs(iRecs,:)) filesep];
        load([cPath 'blueAvg.mat']);

        if exist([cPath 'opts2.mat'],'file')
            load([cPath 'opts2.mat'])
        else
            fprintf('No opts2 found in rec: %s. Trying to align with LEAP model.\n', recs(iRecs,:));
            fprintf('Rec: %i/%i\n',iRecs,size(recs,1));
            opts = rateDisc_alignSingleSession(blueAvg, alignPic, modelPath, alignOpts, cropRange);
            opts.fPath = cPath;
            save([cPath 'opts2.mat'], 'opts');
        end
        
        % keep images and translation values
        rawMovie(:,:,iRecs) = blueAvg;
        try
            snapMovie(:,:,iRecs) = alignAllenTransIm(blueAvg, opts.transParams);
        catch
            opts = alignOpts; %in case LEAP failed badly try initial alignment
            snapMovie(:,:,iRecs) = alignAllenTransIm(blueAvg, opts.transParams);
            save([cPath 'opts2.mat'], 'opts');
        end 
        appShift(:, iRecs) = opts.transParams.tC; %translation
        appRotate(iRecs) = opts.transParams.angleD; %rotation
        appScale(iRecs) = opts.transParams.scaleConst; %scaling
        
    catch
        rejIdx(iRecs) = true;
        fprintf('Failed to align rec: %s\n', recs(iRecs,:));
    end
end

% drop rejected recordings
recs(rejIdx, :) = [];
snapMovie(:,:,rejIdx) = [];
rawMovie(:,:,rejIdx) = [];
appShift(:, rejIdx) = [];
appRotate(rejIdx) = [];
appScale(rejIdx) = [];

%% run alignment tool
rawMovie(isnan(rawMovie)) = 0; %don't use NaNs
BrainAligner(alignPic, rawMovie, appShift, appRotate, appScale); % use this to create appRotate and appShift

%% save down new alignment
alignMovie = NaN(size(alignPic,1), size(alignPic,2), size(recs,1), 'single'); %aligned vessel images
for iRecs = 1 : size(recs,1)
    
    cPath = [fPath strtrim(recs(iRecs,:)) filesep];
    load([cPath 'opts2.mat'], 'opts');
    opts.transParams.tC = appShift(:, iRecs); %translation
    opts.transParams.angleD = appRotate(iRecs); %rotation
    opts.transParams.scaleConst = appScale(iRecs); %scaling
    save([cPath 'opts2.mat'], 'opts');
    
    load([cPath 'blueAvg.mat'], 'blueAvg');
    alignMovie(:,:,iRecs) = alignAllenTransIm(blueAvg, opts.transParams);

end