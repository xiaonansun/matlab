function processCaImagingMCPnev(paramsFileName)
% processCaImagingMCPnev(paramsFileName)
%
% This is the master function for Matt's Ca imaging pathway. It performs
% the following steps:
% 1) Load a parameters file from a folder in the current directory named
%    improcparams/paramsFileName. This file should be written by
%    writeCaProcessParams()
% 2) Motion correction using dftregistration. This is a fast, subpixel,
%    rigid-frame translation. The core function is preprocessCaMovies(). If
%    more than one channel is recorded, motion correction can be done on a
%    particular channel and the results can be applied to other channels.
%    The motion correction step will be skipped automatically if it has
%    already been performed.
%    This step will write tif files ending in _MCM, and mat files named
%    "yymmdd_xxx_zzz", where yymmdd is the date that the mdf file was
%    recorded, xxx and zzz are the tif 'major' and 'minor' numbers,
%    respectively; these mat files include DFToutputs, ie the output of the
%    dftregistration, for each tif 'minor' file. For the entire movie (ie
%    all tif files corresponding to a mdf file) a mat file (named
%    "yymmdd_xxx") will be created which includes badFrames, pmtOffFrames,
%    DFTouputs, and maskBounds corresponding to the entire movie.
% 3) If requested, append to the "yymmdd_xxx" mat file the median, std dev,
%    max, and range images of the entire movie.
% 4) If requested, run Etychios and Liam's automatic source separation
%    algorithm. The core function is applyPnevPaninskiCaSourceSep(). The
%    results will be saved to a mat file named
%    [yymmdd_xxx-PnevPanResults-"nowStr"], where nowStr is the date and
%    time that the algorithm was run.
% 5) If requested, merge the results into the behavior file specified by
%    the supplied params file. This can optionally call a cleanup function
%    (from the params file) on the behavior data first. Saves yet another
%    file.
%
% To evaluate the output, load the result (located in the dataset's folder)
% and see ROIContoursPnev()
%
% Written by Matt Kaufman and Farzaneh Najafi


MCMSuffix = '_MCM';
% maxMaskWidth = 20;


%% Load parameters

fprintf('loading file %s\n', fullfile(pwd, 'improcparams', paramsFileName))
loadVar = load(fullfile(pwd, 'improcparams', paramsFileName));
params = loadVar.params;


%{
%% Rename tif files so tif minor is XXX instead of XX.
% set current tif file names
if any(params.oldTifName)
    files = dir(params.tifFold);
    files = files(~[files.isdir]);

    % Parse all the filenames, save results into files struct array
    for f = 1:length(files)
      [nums, valid] = parseCaImagingTifName(files(f).name);
      files(f).nums = nums;
      files(f).valid = valid;
    end
    files = files([files.valid] & arrayfun(@(f) ismember(f.nums(2), params.tifNums(1,2)), files'));
    
    % rename tif files to the newest mview standard
    for ifile = 1: size(params.tifNums, 1)
        if params.oldTifName(ifile)
            fprintf('Renameing file %s\n', files(ifile).name)
            oldName = files(ifile).name;
            src = fullfile(params.tifFold, oldName);
            
            newName = assembleCaImagingTifName(params.tifNums(ifile,:), params.oldTifName(ifile));
            dest =  fullfile(params.tifFold, newName);
            
            movefile(src, dest)
        end
    end
end
%}



%% Set maxNumCompThreads

warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
% warning('off', 'MATLAB:nargchk:deprecated')
warning('on','setImagedChannels:noImgDesc')

% FN: Choose limit_threads similar to -pe threads in the cluster-submitted job, ie ~16.
fprintf('Original maxNumCompThreads = %.1d\n', maxNumCompThreads)
if params.limit_threads
    origThreadCount0 = maxNumCompThreads;
    maxNumCompThreads(params.limit_threads);
    fprintf('maxNumCompThreads after initial reset = %.1d\n', maxNumCompThreads)
end


%% Figure out if preprocessing has been done for every tif file

preprocDone = true;
if ~params.motionCorrDone
    preprocDone = false;
end


% FN: the part below has issues and needs work. Even when MCM tif exists,
% the raw tif will make params.tifNums(f, 4) equal 0 and hence motion
% correction will be again performed! It is better and simpler to give
% motionCorrDone as an input to writeCaProcessParams, instead of
% having the script set it.
%{
nFiles = size(params.tifNums, 1);
for f = 1:nFiles
    if params.tifNums(f, 4) == 0
        fprintf('At least one file does not have a motion-corrected version. Performing motion correction.\n');
        preprocDone = false;
        break;
    end
    badFramesName = fullfile(params.tifFold, assembleCaImagingTifName(params.tifNums(f, :), 1));
    if ~exist(badFramesName, 'file')
        fprintf('At least one file does not have a _badFrames file. Performing motion correction.\n');
        preprocDone = false;
        break;
    end
end
%}


%% Do preprocessing if it hasn't been done already

% date_major = sprintf('%06d_%03d', params.tifNums(1, 1:2));
% Codes below allow for names like "151102_001_002" in case multiple sessions (ie mdf files, aka tif majors) are being analyzed at the same time.
u = unique(params.tifNums(:,2));
r = repmat('%03d-', 1, length(u)); 
r(end) = [];
date_major = sprintf(['%06d_', r], params.tifNums(1,1), u'); 


if preprocDone
    fprintf('Motion correction has already been performed, skipping\n');
    t1 = tic;
    
    %% Load tifs and set movieMC
    if isempty(params.channelsToRead) % read all saved channels.
        chAll = unique(params.tifNums(:,4)); % channels saved.
        chAll(isnan(chAll)) = [];
    else
        chAll = params.channelsToRead;
        chAll = chAll(:);
    end
    
    movieMC = cell(1, max([chAll', params.dftRegCh, params.activityCh]));
    for ch = chAll'
        % Get list of MCM tif files corresponding to channel ch.
        tifNumsCh = params.tifNums(params.tifNums(:,4)==ch,:);
        tifList = cell(1, size(tifNumsCh,1));
        for f = 1:length(tifList)
            tifList{f} = fullfile(params.tifFold, assembleCaImagingTifName(tifNumsCh(f, :), params.oldTifName(f)));
        end
        
        %       movieMC{ch} = readTiffSet(tifList);
        for t = 1:length(tifList)
            fprintf('Reading tif file %s\n', tifList{t})
            movieMC{ch} = cat(3, movieMC{ch}, bigread2(tifList{t}));
        end
    end
    
    fprintf('\nLoading data took %0.1f s\n\n', toc(t1));
    
    
    %% Set badFrames and pmtOffFrames
    if params.allTifMinors % analyze all tif minors.
        % Load badFrames and pmtOffFrames
        load(fullfile(params.tifFold, date_major), 'badFrames', 'pmtOffFrames')
        
        if ~exist('pmtOffFrames', 'var') % a = whos('-file', '151021_001.mat'); ~any(strcmp('pmtOffFrames', {a.name}))
            pmtOffFrames = cell(size(badFrames));
            for ch = 1:length(pmtOffFrames)
                pmtOffFrames{ch} = false(length(badFrames{ch}),1);
            end
        end
        
        
    else % only analyze some tifMinors.
        % load 'badFramesTif', 'pmtOffFramesTif' for each tifMinor, and
        % concatenate them to set 'badFrames', 'pmtOffFrames'.
        badFrames = cell(1, max(chAll));
        pmtOffFrames = cell(1, max(chAll));
        
        tifMinor = unique(params.tifNums(:,3))';
        for itm = 1:length(tifMinor)
            a = dir(fullfile(params.tifFold, [date_major, '_*', num2str(tifMinor(itm)), '.mat']));
            load(fullfile(params.tifFold, a.name), 'badFramesTif', 'pmtOffFramesTif')
            
            if exist('badFramesTif', 'var')
                for ch = 1:length(badFramesTif)
                    badFrames{ch} = [badFrames{ch}; badFramesTif{ch}];
                    pmtOffFrames{ch} = [pmtOffFrames{ch}; pmtOffFramesTif{ch}];
                end
            end
        end
        
        
        if ~exist('badFramesTif', 'var') % for some days badFramesTif and pmtOffFramesTif were not saved, so you need to set them here.
            % load badFrames and pmtOffFrames corresponding to the entire mdf file.
            load(fullfile(params.tifFold, date_major), 'badFrames', 'pmtOffFrames')
            
            if ~exist('pmtOffFrames', 'var') % a = whos('-file', '151021_001.mat'); ~any(strcmp('pmtOffFrames', {a.name}))
                pmtOffFrames = cell(size(badFrames));
                for ch = 1:length(pmtOffFrames)
                    pmtOffFrames{ch} = false(length(badFrames{ch}),1);
                end
            end
            
            % load DFToutputs of the first tif minor, and use its size to figure
            % out how many frames were saved per Tif file (except for the last
            % tif file that includes whatever frame is remained).
            a = dir(fullfile(params.tifFold, [date_major, '_*01.mat']));
            load(fullfile(params.tifFold, a.name), 'DFToutputs')
            nFramesPerMovie_est = size(DFToutputs{find(cellfun(@(x) ~isempty(x), DFToutputs),1)} , 1);
            
            % set the frames corresponding to tifMinor
            cs = [0:nFramesPerMovie_est:size(badFrames{1},1) length(badFrames{1})]; % [0 cumsum(nFramesPerMovie)];
            frs = [];
            for itm = tifMinor
                frames = cs(itm)+1 : cs(itm+1);
                frs = [frs , frames];
            end
            
            % extract those frames from badFrames and pmtOffFrames
            for ch = 1:length(pmtOffFrames)
                badFrames{ch} = badFrames{ch}(frs);
                pmtOffFrames{ch} = pmtOffFrames{ch}(frs);
            end
        end
    end
    
    
    %%
else
    fprintf('Performing motion correction\n');
    t1 = tic;
    
    % Get list of non-MCM tif files.
    tifNumsNoMCM = params.tifNums(isnan(params.tifNums(:,4)),:);
    tifList = cell(1, size(tifNumsNoMCM,1));
    for f = 1:length(tifList)
        tifList{f} = fullfile(params.tifFold, assembleCaImagingTifName(tifNumsNoMCM(f, :), params.oldTifName(f)));  % tifList does not include any _MCM tif files.
    end
    
    regTif = params.regTifFile;
    regFrameNums = params.regFrameNums;
    dftRegChannel = params.dftRegCh;
    channels2write = params.channelsToWrite;
    maxMaskWidth = params.maxMaskWidth;
    analysisDir = params.analysisFold;    
    
    [movieMC, badFrames, pmtOffFrames] = preprocessCaMovies(tifList, regTif, regFrameNums, dftRegChannel, channels2write, MCMSuffix, maxMaskWidth, analysisDir, params.pmtOffThr, params.makeMCMrepMovie, params.frsExclude);
    
    paramsRegist.regTifFile = params.regTifFile;
    paramsRegist.regFrameNums = params.regFrameNums;
    save(fullfile(params.tifFold, date_major), 'paramsRegist', '-append')
    fprintf('\nMotion correction took %0.1f s\n\n', toc(t1));
end

ind = find(cellfun(@(x)~isempty(x), movieMC), 1);
[imHeight, imWidth, nFrames] = size(movieMC{ind});  % size(movieMC{params.activityCh});
save(fullfile(params.tifFold, date_major), 'params', 'imHeight', 'imWidth', '-append') % This mat file is already created during motion correction, which is why we append to it here.


%% Convert movieMC to single 
% OK to convert to double here, because applyPnevPaninskiCaSourceSep checks
% to see whether it's already a double and will skip conversion if it is.
% Also std requires single or double inputs.
if params.pnevActivity || params.saveGoodMovieStats
    t1 = tic;
    for ch = 1:length(movieMC)
        if ~isa(movieMC{ch}, 'single')
            fprintf('Converting to single, channel %d movie...\n', ch);
            tic;
            movieMC{ch} = single(movieMC{ch});
            fprintf('%0.2fs\n', toc);
        end
    end
    fprintf('\nConverting movie to single took %0.1f s\n\n', toc(t1));
end


%% Produce and save averaged images: Median, std dev, max, range

if params.saveGoodMovieStats % && ~preprocDone % in writeCa we make the default value of saveGoodMovieStats 0 if motion correction is done. 
    
    t1 = tic;
    % Pre-allocate cell arrays
    aveImage = cell(1,length(movieMC));
    medImage = cell(1,length(movieMC));
    sdImage = cell(1,length(movieMC));
    maxImage = cell(1,length(movieMC));
    rangeImage = cell(1,length(movieMC));
    
    for ch = 1:length(movieMC)
        if ~isempty(movieMC{ch})
            fprintf('Generating various projection images for channel %d...\n', ch);
            tic;
            
            % Pre-allocate
            aveImage{ch} = NaN(imHeight, imWidth);
            medImage{ch} = NaN(imHeight, imWidth);
            sdImage{ch} = NaN(imHeight, imWidth);
            maxImage{ch} = NaN(imHeight, imWidth);
            rangeImage{ch} = NaN(imHeight, imWidth);
            
            % Compute images pixel by pixel. We do this so that we can skip bad frames
            % without having to allocate another entire movie.
            for i1 = 1:imHeight
                for i2 = 1:imWidth
                    pixel = movieMC{ch}(i1, i2, ~badFrames{ch} & ~pmtOffFrames{ch});
                    
                    aveImage{ch}(i1, i2) = mean(pixel);
                    medImage{ch}(i1, i2) = median(pixel);
                    sdImage{ch}(i1, i2) = std(pixel);
                    maxImage{ch}(i1, i2) = max(pixel);
                    rangeImage{ch}(i1, i2) = range(pixel);
                end
            end
            
            fprintf('%0.2fs\n', toc);
        end
    end
    clear pixel
    
    % Save results
    a = matfile(fullfile(params.tifFold, [date_major, '.mat'])); % check vars saved in date_major mat file.
    if isprop(a, 'medImage') && any(cellfun(@length, a.medImage)) % if average images were already saved for a channel, and now we are getting them for a different channel, make sure we don't overwrite the previous variables.
        if isfield(params, 'channelsToRead')
          channels = params.channelsToRead;
        else
          channels = params.channelsToWrite;
        end
        save(fullfile(params.tifFold, [date_major, '_ch', num2str(channels)]), 'aveImage', 'medImage', 'sdImage', 'maxImage', 'rangeImage')
    else
        save(fullfile(params.tifFold, date_major), 'aveImage', 'medImage', 'sdImage', 'maxImage', 'rangeImage', '-append')
    end
    
    fprintf('\nComputing and saving summary images took %0.1f s\n\n', toc(t1));
end


%%%%%%%%%%%%%%%%
% Don't delete the commented codes below. Uncomment it if you need to
% only compute aveImage (for those days that don't have this variable
% saved.)

if params.onlyAveImg
    aveImage = cell(1,length(movieMC));
    for ch = 1:length(movieMC)
        if ~isempty(movieMC{ch})
            fprintf('Generating various projection images for channel %d...\n', ch);
            tic;
            aveImage{ch} = NaN(imHeight, imWidth);
            for i1 = 1:imHeight
                for i2 = 1:imWidth
                    pixel = movieMC{ch}(i1, i2, ~badFrames{ch} & ~pmtOffFrames{ch});

                    aveImage{ch}(i1, i2) = mean(pixel);

                end
            end

            fprintf('%0.2fs\n', toc);
        end
    end
    clear pixel
    save(fullfile(params.tifFold, date_major), 'aveImage', '-append')    
end


%% Run Eftychios and Liam's algorithm (CNMF) to identify ROIs and set their (deconvolved) fluorescent traces.

if params.pnevActivity
    
    if ~params.cnmfStep1Done % run CNMF from the beginning
        nowStr = datestr(now, 'yymmdd-HHMMSS');
        pnevFileName = fullfile(params.tifFold, [date_major, '_ch', num2str(params.activityCh), '-PnevPanResults-', nowStr]);

        pnev_inputParams.K = params.numComps;                                % number of components to be found
        pnev_inputParams.temp_sub = params.tempSub;                          % temporal subsampling for greedy initiation, set to 1 for no down sampling.
        pnev_inputParams.space_sub = params.spaceSub;                        % spatial subsampling for greedy initiation, set to 1 for no down sampling.
        pnev_inputParams.deconv_method = params.deconv_method;               % activity deconvolution method ('constrained_foopsi', 'MCMC')
    %     pnev_inputParams.finalRoundMCMC = params.finalRoundMCMC;             % do a final round of MCMC method (if false, after merging 2 iterations of const foopsi will be done. If true, after merging 1 iter of const foopsi and 1 iter of MCMC will be done.)
        pnev_inputParams.doPlots = params.doPlots;                           % if true, some figures and a movie will be made.
        %     pnev_inputParams.parallelTempUpdate = params.parallelTempUpdate;     % do parallel temporal updating.
        pnev_inputParams.save4debug = params.save4debug;                     % save Eftychios's variables (eg A,C,etc) after each step for debug purposes.
        pnev_inputParams.MCMC_B = params.MCMC_B;                             % number of burn in samples. eg. 200
        pnev_inputParams.MCMC_Nsamples = params.MCMC_Nsamples;               % number of samples after burn in. eg. 200
        pnev_inputParams.MCMC_prec = params.MCMC_prec;                       % specifies the extent of discarding the long slowly decaying tales of the ca response. eg. 1e-2
        pnev_inputParams.save_merging_vars = params.save_merging_vars;       % whether to save some vars related to merging components.
        pnev_inputParams.search_dist = params.search_dist;                   % search distance when updating spatial components.
        pnev_inputParams.init_method = params.init_method;
        pnev_inputParams.temporal_parallel = params.temporal_parallel;
        pnev_inputParams.limit_threads = params.limit_threads;
        pnev_inputParams.orderROI_extractDf = params.orderROI_extractDf;
        pnev_inputParams.p = params.ARmodelOrder;                          % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
        pnev_inputParams.multiTrs = params.multiTrs;
        pnev_inputParams.maxFrsForMinPsn = params.maxFrsForMinPsn;         % maximum number of frames of the movies for computing min(Y) and P.sn. It helps with speed but assumes that the 1st maxFrsForMinPsn frames give a good estimate for min(Y) and P.sn, so no need to look at the entire movie.
        pnev_inputParams.poolsize = params.poolsize;                       % number of parallel workers. If set to 0, the default NumWorkers of the cluster Pool will be used. 
        pnev_inputParams.doMerging = params.doMerging;                     % whether to merge components (which will be followed by another round of spatial and temporal updating).
        pnev_inputParams.noise_norm = params.noise_norm;                   % whether to do normalization by noise estimate prior to initialization. (default: 0) 
        pnev_inputParams.noise_norm_prctile = params.noise_norm_prctile;   % minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
        pnev_inputParams.bas_nonneg = params.bas_nonneg;                   % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
        pnev_inputParams.customMerge = params.customMerge;                 % flag for using Farzaneh's function mergeROIs_set to define merged_ROIs.
        pnev_inputParams.cnmfStep1Done = params.cnmfStep1Done;             % if 1, only deconvolution is remained, otherwise run CNMF from the beginning.

        pnev_inputParams.tau = params.tau; %4;                              % std of gaussian kernel (size of neuron)
        pnev_inputParams.merge_thr = 0.8;                                   % merging threshold    
        pnev_inputParams.temp_iter = 2;                                     % number of block-coordinate descent steps
        pnev_inputParams.fudge_factor = .98;                                % bias correction for AR coefficients
    
    else % CNMF is already run, only do the deconvolution step.
        pnevFileName = params.pnevFileName;
        load(pnevFileName) % load all the CNMF variables that are already saved.
%         pnev_inputParams.pnevFileName = params.pnevFileName;
    end
    
    
    %% Run CNMF
    
    if ~params.cnmfStep1Done
        
        if isfield(params, 'brightenNorming') && params.brightenNorming
            if ~exist('medImage', 'var')
                %             clear medImage   % unnecessary, we just checked that medImage isn't defined
                a = matfile(fullfile(params.tifFold, [date_major, '.mat'])); % check vars saved in date_major mat file.
                if isprop(a, 'medImage')
                    load(fullfile(params.tifFold, date_major), 'medImage')
                else
                    medImage{params.activityCh} = median(movieMC{params.activityCh}, 3);
                end
            end
            normingMedianImage = brightenFilter2DGauss(medImage{params.activityCh},params);
            pnev_inputParams.normingMedianImage = normingMedianImage;
        end

        [A, C, S, C_df, Df, b, f, srt, srt_val, Ain, options, P, merging_vars, YrA, Yr, Cin, bin, fin, nA,AY,AA] = demo_script_modif(movieMC{params.activityCh}, pnev_inputParams);

    end
    
    
    %% Do deconvolution for the multi-trial case (FN)
    
    if params.multiTrs 
        
        load(fullfile(params.tifFold, date_major), 'cs_frtrs', 'Nnan_nanBeg_nanEnd')
        if exist('cs_frtrs', 'var')
            params.cs_frtrs = cs_frtrs;
        end 
        if exist('Nnan_nanBeg_nanEnd', 'var')
            params.Nnan_nanBeg_nanEnd = Nnan_nanBeg_nanEnd;
        end            

        P.p = 2;
        [A, C, S, C_df, Df, f, P, srt, srt_val, nA, YrA, AY] = update_tempcomps_multitrs(C, f, A, b, AY, AA, P, options, params); %, Yr) %(C, f, A, b, Yr, P, options, params);
        
%         save(fullfile(params.tifFold, [date_major, '_ch', num2str(params.activityCh), '-PnevPanResults-', nowStr]), ...
%             'A', 'C', 'S', 'C_df', 'Df', 'f', 'P', 'srt', 'srt_val', 'nA', 'YrA', 'AY', '-append')
    end    
%     P = rmfield(P, 'psdx');  % FN commented bc cluster_pixels is false by default in the current version. % MTK: psdx is huge and not very useful %

    
    %% Save the results
    
    if ~params.cnmfStep1Done
        fprintf('Saving Pnev results.\n')        
        save(pnevFileName, ... % we need version 7.3 or it will skip saving some vars.
            'A', 'C', 'S', 'C_df', 'Df', 'b', 'f', 'srt', 'srt_val', 'Ain', 'options', 'P', 'pnev_inputParams', 'merging_vars', 'bin', 'fin', 'YrA', 'nA', 'AY', 'AA', '-v7.3');
    
    else % Rewrite the new variables 
        fprintf('Appending Pnev results.\n')        
        save(pnevFileName, '-append', ...
            'A', 'C', 'S', 'C_df', 'Df', 'f', 'P', 'srt', 'srt_val', 'nA', 'YrA', 'AY');        
    end
    

end


%% Compute manual DF/F for CNMF-found ROIs

% Set pnevFileName if it is not in the workspace
if ~exist('pnevFileName', 'var')
    nameTemp = [date_major, '_ch', num2str(params.activityCh),'-Pnev*'];
    nameTemp = dir(fullfile(params.tifFold, nameTemp));

    if ~isempty(nameTemp) % pnevFileName is already saved    
        [~,i] = sort([nameTemp.datenum], 'descend'); % in case there are a few pnev files, choose which one you want!
        % disp({pnevFileName(i).name}')
        pnev2load = 1; % use the most recent file.
        pnevFileName = nameTemp(i(pnev2load)).name;
        pnevFileName = fullfile(params.tifFold, pnevFileName);
    end
end


if params.manualActivity % FN: whether to compute manual activity (ie average pixel intensity over time) for ROIs found by Eftychios's algorithm.
    if ~exist('A', 'var')
        load(pnevFileName, 'A')
    end
    for ch = 1:length(movieMC)
        if ~isempty(movieMC{ch})
            fprintf('Computing manual activity for CNMF ROIs, channel %d', ch)
            tic;
            activity_man_eftMask = manualROIactivityFromEftMaskMovieMC(movieMC{ch}, A, imHeight, imWidth);

            fprintf('... took %.1f s\n', toc)

            nam2save = ['activity_man_eftMask_ch', num2str(ch)]; % activity_man_eftMask_ch2
            eval([nam2save ' = activity_man_eftMask;'])
            save(pnevFileName, ...
                '-append', nam2save) %, 'activity_man_eftMask')
        end
    end
end    
%     activity = C;    
    %{
    [spatialBasis, activity, spatialBackground, backgroundActivity, spiking, Df, activityDf, spikingDf, pnevParam, ROIOrderMeasure, greedySpatialBasis] = ...
        applyPnevPaninskiCaSourceSep(movieMC{params.activityCh});
%     [A, C, b, f, S, Df, C_df, S_df, P, ROIOrderMeasure, Ain]

    save(fullfile(params.tifFold, [date_major, '-PnevPanResults-', nowStr]), ...
        'spatialBasis', 'activity', 'spatialBackground', 'backgroundActivity', ...
        'spiking', 'Df', 'activityDf', 'spikingDf', 'pnevParam', 'ROIOrderMeasure', 'greedySpatialBasis');
    %}


%% Compute correlation between A and raw hightlight-reel of each ROI

% Our version
if params.performCorr
    if ~exist('S', 'var')
        load(pnevFileName, 'A', 'S')
    end
    [highlightCorrROI,highlightPatchAvg,roiPatch,highlightTifs] = CorrelationThresholding(A,S,movieMC{params.activityCh},params.padCorrelation);

    if params.saveHighlight
        save(pnevFileName, ...
            '-append', 'highlightCorrROI','highlightPatchAvg','roiPatch', 'highlightTifs');
    else
        save(pnevFileName, ...
            '-append', 'highlightCorrROI','highlightPatchAvg','roiPatch');
    end

end


%%%% Efty's version
if params.performCorrEP
    if ~exist('C', 'var')
        load(pnevFileName, 'A', 'C', 'b', 'f')
    end

    [d1,d2,T] = size(movieMC{params.activityCh});                                % dimensions of dataset
    d = d1*d2;                                          % total number of pixels
%     Yr = reshape(movieMC{params.activityCh} - min(movieMC{params.activityCh}(:)),d,T);

    %%%
    options.space_thresh = .4;
    options.time_thresh = .4;
    options.A_thresh = .1;
    options.Npeaks = 50; % 10;
    options.peak_int = 0:1;% -5:24;

    %%%
%     [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(Yr,A,C,b,f,options);
    [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(reshape(movieMC{params.activityCh} - min(movieMC{params.activityCh}(:)),d,T),A,C,b,f,options);

    %%%
    save(pnevFileName, ...
        '-append', 'rval_space','rval_time','ind_space','ind_time');        
end



%% Reset maxNumCompThreads to its original value

if params.limit_threads
    maxNumCompThreads(origThreadCount0);
end
fprintf('maxNumCompThreads, final reset = %.1d\n', maxNumCompThreads) % FN check for num threads.






%% If requested, merge into behavior

if ~isempty(params.behavFile)
    
    %% Get trialization info
    
    framesPerTrial = cell(1, length(params.binFiles));
    trialNumbers = cell(1, length(params.binFiles));
    frame1RelToStartOff = cell(1, length(params.binFiles));
    
    for f = 1:length(params.binFiles)
        [framesPerTrial{f}, trialNumbers{f}, frame1RelToStartOff{f}] = ...
            framesPerTrialStopStart3An(params.binFiles{f}, params.framecountFiles{f}, params.headerBug);
    end
    
    framesPerTrial = [framesPerTrial{:}];
    trialNumbers = [trialNumbers{:}];
    frame1RelToStartOff = [frame1RelToStartOff{:}];
    
    
    %% Load the alldata file
    
    loadVar = load(params.behavFile);
    alldata = loadVar.all_data;
    
    
    %% If specified, call cleanup function
    
    if ~isempty(params.behavFcn)
        alldata = feval(params.behavFcn, alldata);
    end
    
    
    %% Merge activity into alldata, save
    
    ad = mergeActivityIntoAlldataPnev(alldata, C_df(1:end-1, :)', S_df', framesPerTrial, trialNumbers, frame1RelToStartOff, badFrames{params.activityCh}, pmtOffFrames{params.activityCh});
    %   ad = mergeActivityIntoAlldataPnev(alldata, activityDf', spikingDf', framesPerTrial, trialNumbers, frame1RelToStartOff, badFrames{params.activityCh}, pmtOffFrames{params.activityCh});
    save(params.mergedName, 'ad');
    
end


fprintf('Done.\n');
