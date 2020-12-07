function params = writeCaProcessParams(outName, mouse, imagingFold, tifMajor, P)
% params = writeCaProcessParams(outName, mouse, imagingFold, tifMajor, P)
%
% Goal: make it easy to write a parameters file for use with
% processCaImagingMCPnev()
% 
% INPUTS
%
% outName       -- string, name of parameters file to output. Will be
%                  placed in ./improcparams/
% mouse         -- string, name of mouse
% imagingFold   -- string, name of folder containing imaging data relative
%                  to dataPath/imaging
% tifMajor      -- tif "major" numbers, from MScan. May have length > 1
%
% P             -- structure with the following fields. All are optional,
%                  except for signalCh.
% signalCh      -- channel whose signal activity we want to analyse (usually gcamp channel)
% regFileNums   -- optional, 1 x 3 numeric array, containing info about which file to
%                  use for motion correction. Should be: [major minor
%                  channel]. Registration will be skipped if empty or NaN.
% regFrameNums  -- optional, numeric array or cell array. If numeric array, it
%                  indicates frames to use from reg file for motion
%                  correction. If cell array it indiates trial(s) to use
%                  for motion correciton. For simplicity, this trial must
%                  be chosen from the 1st tif file, and regTif must be the
%                  1st tif file (ie minor=1).
% behavName     -- optional. Name of behavioral file to merge into. Ignored
%                  if empty string (i.e., merging will not be attempted if
%                  not specified or empty)
% behavFcn      -- optional. Function to be called on behavioral data for
%                  cleanup and merge. Ignored if empty string.
% headerBug     -- optional, default 0. Whether the analog channel was from
%                  the buggy version of MScan that messes up the header.
% maxMaskWidth  -- optional, default 20. Max number of pixels that will be
%                  masked from each side of a frame after motion
%                  correction. Frames that required more pixel shifts in
%                  order to be registered will be marked as badFrames. Set
%                  to inf to allow masking as much as needed by pixelShifts.
% channelsToWrite -- optional. Movie from these channels will be written to
%                  tif files. If not provided or if empty, both signalCh
%                  and dftRegCh channels will be written.
% analysisFolder  -- optional. If 1, files resulting from the analysis will
%                  be saved in a separate folder named "analysis". If not
%                  provided or 0, they will be saved in the same folder as
%                  the
%                  imaging folder.
% motionCorrDone -- optional. If 1, indicates motion correction files are
%                  already saved and it wont be performed again.
% pmt_th         -- optional. If provided, it determins the average pixel
%                  threshold for finding pmtOff frames (ie frames during
%                  which PMT was off). The pixel values of these frame will
%                  be turned to NaN (since movieMC is uint16, NaN will be
%                  converted to 0, so pmtOff frames will be all 0). If not
%                  provided, no frame will be marked as pmtOffFrame. 
%                  If is a cell it indicates the index of pmtOffFrames.
% channelsToRead -- Channels to load in case of motionCorrDone. Default:
%                  all saved channels.
% saveGoodMovieStats -- Whether or not to save stats (max, sd, median,
%                  range) of goodMovie (movie ignoring badFrames and
%                  pmtOffFrames). Default: true if motion correction not performed, otherwise false.
% pnevActivity -- optional, default true, indicating that Eftychios's
%                  algorithm will be used to drive ativity traces. If
%                  provided and 0, Eftychios's algorithm will not be executed.
% tifMinor     -- optional, default all tif files of a mdf file will be
%                  included in params. If provided and non-empty, only
%                  tifMinor files will be included in params for analysis.
% saveParams   -- optional, default true. If true, params will be saved to
%                  improcess folder.
%
% The following optional fields are related to Eftychios's algorithm
% numComps           -- number of components to be found, default:200
% tempSub            -- temporal subsampling for greedy initiation, set to 1 for no down sampling. default:3
% spaceSub           -- spatial subsampling for greedy initiation, set to 1 for no down sampling. default:2
% deconv_method      -- default: 'constrained_foopsi'; activity deconvolution method ('constrained_foopsi', 'MCMC', ...). If MCMC, after merging 1 round of const foopsi and then 1 round of MCMC will be done.
% doPlots            -- make some figures and a movie, default:false
% save4debug         -- save Eftychios's variables (eg A,C,etc) after each step for debug purposes. default:false
% MCMC_B             -- MCMC deconvolution method, number of burn in samples. default 300
% MCMC_Nsamples      -- MCMC deconvolution method, number of samples after burn in. default 400
% MCMC_prec          -- MCMC deconvolution method, the extent of discarding the long slowly decaying tales of the ca response. default 1e-2
% save_merging_vars  -- whether to save some vars related to merging components; allows offline assessment if desired.
% search_dist        -- search distance when updating spatial components.
% init_method        -- whether to use 'greedy' or 'sparse_NMF' to initialize spatial components
% limit_threads      -- whether to limit the number of threads in
%                       processCaImagingPCPnev and demo_script_modif.
%                       limit_threads should be close to -pe threads when
%                       submitting jobs to the cluster, eg 16. If >0,
%                       specifies maxNumCompThreads, ie *total* number of
%                       threads. For parfor loops, limit_threads will be
%                       divided among numWorkers (specified by poolsize),
%                       so if limit_threads=16, and poolsize=8, outside
%                       parfor maxNumCompThreads will be 16, and inside
%                       parfor it will be 16/8=2 threads per worker.)
% temporal_parallel  -- whether to parallelize temporal updating
% maxFrsForMinPsn    -- default: []; maximum number of frames of the movies
    % for computing min(Y) and P.sn. It helps with speed but assumes that the
    % 1st maxFrsForMinPsn frames give a good estimate for min(Y) and P.sn, so
    % no need to look at the entire movie. by default all frames will be used.
% poolsize           -- number of parallel workers. If set to 0, the default NumWorkers of the cluster Pool will be used. 
% manualActivity     -- (default 0); whether to compute manual activity (ie average pixel intensity over time) for ROIs found by Eftychios's algorithm.
% doMerging          -- (default 1); whether to use Eftychios's merge_components (which will be followed by another round of spatial and temporal updating).
% noise_norm         -- (default 0) whether to do normalization by noise estimate prior to initialization. 
% noise_norm_prctile -- (default 2) minimum noise level (as percentile of P.sn) used in the normalization prior to initialization 
% bas_nonneg         -- (default 1) flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) 
% customMerge        -- (default 0) if 1, Farzaneh's function mergeROIs_set will be used to define merged_ROIs.
% makeMCMrepMovie    -- (default 0) if 1, write a movie of random frames (after motion correction) to assess how good the motion correction was. 
%
% Written by Matt Kaufman and Farzaneh Najafi


%% Feedback, setup

fprintf('Mouse: %s\n', mouse)

if isunix
    if isempty(strfind(pwd, 'grid')) %isempty(strfind(pwd, 'sonas')) % Farzaneh's Unix in the office
        dataPath = '~/Shares/Churchland_nlsas_data/data'; '~/Shares/Churchland/data';
        altDataPath = '~/Shares/Churchland_hpc_home/space_managed_data';
    else
        dataPath = '/sonas-hs/churchland/nlsas/data/data';
        altDataPath = '/sonas-hs/churchland/hpc/home/space_managed_data';
    end
elseif ispc
    dataPath = '\\sonas-hs.cshl.edu\churchland\data'; % FN
end


%% Imaging tif files

params.tifFold = fullfile(dataPath, mouse, 'imaging', imagingFold);

if ~exist(params.tifFold, 'dir') && exist('altDataPath', 'var')
  params.tifFold = fullfile(altDataPath, mouse, 'imaging', imagingFold);
end


if ~isfield(P, 'analysisFolder') || ~P.analysisFolder
    params.analysisFold = params.tifFold;
elseif P.analysisFolder
    params.analysisFold = fullfile(dataPath, mouse, 'analysis', imagingFold);
end

%{
imfilename = sprintf('%s_%03d_*.TIF*', imagingFold, tifMajor);
files = dir(fullfile(params.tifFold, imfilename));
files = files(cellfun(@(x)ismember(length(x),[17,25]) && ~isnan(str2double(x(12:13))), {files.name})); % make sure tif file name is of format : YYMMDD_mmm_nn.TIF or YYMMDD_mmm_nn_ch#_MCM.TIF
%}

files = dir(params.tifFold);
files = files(~[files.isdir]);
renameTifAddMinorAll = zeros(size(files)); % FN: flag for renaming tif files that lack a tif minor (happens when an mdf file has <4089 frames).
% Parse all the filenames, save results into files struct array
for f = 1:length(files)
    [nums, valid, renameTif, renameTifAddMinor] = parseCaImagingTifName(files(f).name);
    files(f).nums = nums;
    files(f).valid = valid;
    files(f).renameTif = renameTif;
    renameTifAddMinorAll(f) = renameTifAddMinor;
end

% Take care of renaming tif files without a tif minor (FN).
for i = find(renameTifAddMinorAll')
    o = fullfile(params.tifFold, files(i).name);
    [~,fn,fe] = fileparts(files(i).name);
    n = fullfile(params.tifFold, [fn,'_001',fe]);
    movefile(o, n);
end
   

% Find the valid files, with the right major number
files = files([files.valid] & arrayfun(@(f) ismember(f.nums(2), tifMajor), files'));

if isfield(P, 'tifMinor') && ~isempty(P.tifMinor)
    files = files(arrayfun(@(f) ismember(f.nums(3), P.tifMinor), files'));
    params.allTifMinors = false; % all tif minors are included in params.
else
    params.allTifMinors = true;
end

if isempty(files)
    error('No files found!');
end

% Extract the numbers from each tif filename
params.tifNums = NaN(length(files), 4);
for f = 1:length(files)
    params.tifNums(f, :) = parseCaImagingTifName(files(f).name); % [date, major, minor, channel]
end

params.oldTifName = [files.renameTif];

% Display
fprintf('Found %d tif files. Showing filenames:\n', length(files));
for f = 1:length(files)
    fprintf('%s\n', files(f).name);
end

fprintf('Looking for signal on channel %d\n', P.signalCh);

params.activityCh = P.signalCh;


%% Registration-related info

if ~isfield(P, 'motionCorrDone') || ~P.motionCorrDone
    motionCorrDone = false;
elseif P.motionCorrDone
    motionCorrDone = true;
end
params.motionCorrDone = motionCorrDone;

if ~isfield(P, 'maxMaskWidth')
    params.maxMaskWidth = 20;
else
    params.maxMaskWidth = P.maxMaskWidth;
end


if motionCorrDone || ~isfield(P,'regFileNums') || isempty(P.regFileNums) || any(isnan(P.regFileNums))
    if ~isfield(P, 'channelsToRead') % Channels to load in case of motionCorrDone. Default: all saved channels
        params.channelsToRead = [];
    else
        params.channelsToRead = P.channelsToRead;
    end
    
    params.dftRegCh = P.signalCh; % channel to perform dftregistration on.
    fprintf('Registration info not specified. Assuming motion correction already completed. If not, will cause error only when running main processing\n');
    
else
    
    params.dftRegCh = P.regFileNums(3); % channel to perform dftregistration on.
    
    % FN: if the trial number without motion is provided (instead of the frames
    % without motion), find its corresponding frames. (you need to add the
    % following directory to your path for this to work: repoland\utils\Farzaneh)
    if iscell(P.regFrameNums)
        noMotionTrs = P.regFrameNums{1};
        file2read = fullfile(dataPath, mouse, 'imaging', imagingFold, sprintf('framecounts_%03d.txt', P.regFileNums(1)));
        P.regFrameNums = frameNumsSet(file2read, noMotionTrs);
        %     regFrameNums(regFrameNums > length(tifInfo)/length(channelsSaved)) = []; % regTif must be the 1st tif file.
    end
    
    
    if length(P.regFileNums) ~= 3
        error('regFileNums should have length 3 (tif major number, minor number, and channel)');
    end
    
    regTifFile = assembleCaImagingTifName([params.tifNums(1,1) P.regFileNums(1:2) NaN], params.oldTifName(P.regFileNums(2))); % regFileNums(2) is the tif minor name.
    fprintf('Using channel %d of %s for motion registration\n', P.regFileNums(3), regTifFile);
    
    params.regTifFile = fullfile(params.tifFold, regTifFile);
    %   fullfile(tifFold, sprintf('%s_%03d_%02d_ch%d.TIF', tifStem, regFileNums(1), regFileNums(2), regFileNums(3)));
    params.regFrameNums = P.regFrameNums;
    
    % Check file exists
    if ~exist(params.regTifFile, 'file')
        error('Tif file for registration does not exist: %s', params.regTifFile);
    end
    
    % channelsToWrite
    if ~isfield(P, 'channelsToWrite') || isempty(P.channelsToWrite)
        params.channelsToWrite = unique([params.dftRegCh, P.signalCh]);
    else
        params.channelsToWrite = P.channelsToWrite;
    end
    
    % threshold for pmtOffFrames
    if ~isfield(P, 'pmt_th') || isempty(P.pmt_th)
        params.pmtOffThr = NaN;
    else
        params.pmtOffThr = P.pmt_th; % if the average pixel intensity of a frame is below pmt_th that frame will be marked as a pmtOffFrame and its pixel values will be turned to NaN (since movieMC is uint16, NaN will be converted to 0, so pmtOffFrames will be all 0).
    end
    
end


%%
if ~isfield(P, 'pnevActivity')
    params.pnevActivity = true;
else
    params.pnevActivity = P.pnevActivity;
end


%%
if ~isfield(P, 'saveGoodMovieStats') % Whether to save stats (max, sd, median, range) of goodMovie (movie ignoring badFrames and pmtOffFrames or not). Default: true if motion correction not performed, otherwise false.
    if params.motionCorrDone
        params.saveGoodMovieStats = false;
    else
        params.saveGoodMovieStats = true;
    end
else
    params.saveGoodMovieStats = logical(P.saveGoodMovieStats);
end


%% Merged filename

mergedName = [mouse '_' imagingFold '_Eft'];
params.mergedName = fullfile(dataPath, mouse, 'merged', mergedName);


%% headerBug

if isfield(P, 'headerBug')
    params.headerBug = P.headerBug;
else
    params.headerBug = 0;
end


%% ========  Behavior-related portion  ==========

params.binFiles = {};
params.framecountFiles = {};
params.behavFile = '';
params.behavFcn = '';

if isfield(P, 'behavName') && ~isempty(P.behavName)
    
    %% Behavior file for merging
    
    params.behavFile = fullfile(dataPath, mouse, 'behavior', P.behavName);
    
    
    %% Analog channel files, framecounts files
    
    for m = tifMajor
        params.binFiles{end+1} = fullfile(params.tifFold, sprintf('%d_%03d.bin', params.tifNums(1, 1), m));
        params.framecountFiles{end+1} = fullfile(params.tifFold, sprintf('framecounts_%03d.txt', m));
    end
    
    
    %% Verify that files exist
    
    if ~exist(params.behavFile, 'file') && ~exist([params.behavFile '.mat'], 'var')
        error('Behavior file does not exist: %s', params.behavFile);
    end
    
    for f = 1:length(params.binFiles)
        if ~exist(params.binFiles{f}, 'file')
            error('Analog channel binary file does not exist: %s', params.binFiles{f});
        end
        if ~exist(params.framecountFiles{f}, 'file')
            error('Framecounts file does not exist: %s', params.framecountFiles{f});
        end
    end
    
end


if isfield(P, 'behavFcn') && ~isempty(P.behavFcn)
    params.behavFcn = P.behavFcn;
end


%% Eftychio's algorithm for identifying ROIs and activity
if ~isfield(P, 'numComps')
    P.numComps = 200;
end
params.numComps = P.numComps;


if ~isfield(P, 'tempSub')
    P.tempSub = 3;
end
params.tempSub = P.tempSub;


if ~isfield(P, 'spaceSub')
    P.spaceSub = 2;
end
params.spaceSub = P.spaceSub;


if ~isfield(P, 'temporal_parallel')
    P.temporal_parallel = ~isempty(which('parpool')); 
end
params.temporal_parallel = P.temporal_parallel;


if ~isfield(P, 'deconv_method')
    P.deconv_method = 'constrained_foopsi';
end
params.deconv_method = P.deconv_method;


if ~isfield(P, 'doPlots')
    P.doPlots = false;
end
params.doPlots = P.doPlots;


if ~isfield(P, 'save4debug')
    P.save4debug = false;
end
params.save4debug = P.save4debug;


if ~isfield(P, 'MCMC_B')
    P.MCMC_B = 300;                             
end
params.MCMC_B = P.MCMC_B;


if ~isfield(P, 'MCMC_Nsamples')
    P.MCMC_Nsamples = 400;               
end
params.MCMC_Nsamples = P.MCMC_Nsamples;


if ~isfield(P, 'MCMC_prec')
    P.MCMC_prec = 1e-2;                       
end
params.MCMC_prec = P.MCMC_prec;


if ~isfield(P, 'save_merging_vars')
    P.save_merging_vars = false;
end
params.save_merging_vars = P.save_merging_vars;


if ~isfield(P, 'search_dist')
    P.search_dist = 3;
end
params.search_dist = P.search_dist;


if ~isfield(P, 'limit_threads') % Threading
  P.limit_threads = 0;
end
params.limit_threads = P.limit_threads;


if ~isfield(P, 'init_method')
  P.init_method = 'greedy';
end
params.init_method = P.init_method;


if ~isfield(P, 'orderROI_extractDf')
    P.orderROI_extractDf = true;
end
params.orderROI_extractDf = P.orderROI_extractDf;


if ~isfield(P, 'ARmodelOrder')
    P.ARmodelOrder = 2;  % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
end
params.ARmodelOrder = P.ARmodelOrder;


if ~isfield(P, 'multiTrs')
    P.multiTrs = false;
end
params.multiTrs = P.multiTrs;

%
if params.pnevActivity && params.multiTrs
    for ise = 1:length(tifMajor)
        aa = sprintf('%d_%03d.mat', params.tifNums(1, 1), tifMajor(ise));
        a = fullfile(params.tifFold, aa);
        if exist(a, 'file')==2
            fs = whos('-file', a);
            if ~ismember({'framesPerTrial'}, {fs.name}) % sum(ismember({fs.name}, {'Nnan', 'cs_frtrs'}))~=2
                error('Save framesPerTrial in %s to run Eftys algorithm in the multi-trial mode!', aa)
                %         error('You need to save Nnan and cs_frtrs (in the "date_major" .mat file) to run Eftys algorithm in the multi-trial mode!')
            end
        else
            error('%s does not exist!', aa)
        end
    end
end
%}

if ~isfield(P, 'maxFrsForMinPsn')
    P.maxFrsForMinPsn = [];  % min(Y) and P.sn will be computed on the entire movie.       
end
params.maxFrsForMinPsn = P.maxFrsForMinPsn; % min(Y) and P.sn will be computed on the first maxFrsForMinPsn frames.


if ~isfield(P, 'poolsize') % number of parallel workers. If set to 0, the default NumWorkers of the cluster Pool will be used. 
    P.poolsize = 0;
end
params.poolsize = P.poolsize;


if ~isfield(P, 'manualActivity')
    P.manualActivity = 0;
end
params.manualActivity = P.manualActivity;


if ~isfield(P, 'doMerging')
    P.doMerging = 1; % (default 1); whether to use Eftychios's merge_components (which will be followed by another round of spatial and temporal updating).
end
params.doMerging = P.doMerging;


if ~isfield(P, 'noise_norm') 
    P.noise_norm = 0; % whether to do normalization by noise estimate prior to initialization. (default: false) 
end
params.noise_norm = P.noise_norm;


if ~isfield(P, 'noise_norm_prctile') % it will be only effective if noise_norm is 1.
    P.noise_norm_prctile = 2; % minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
end
params.noise_norm_prctile = P.noise_norm_prctile;


if ~isfield(P, 'bas_nonneg')
    P.bas_nonneg = 1; % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
end
params.bas_nonneg = P.bas_nonneg;


if ~isfield(P, 'customMerge')
    P.customMerge = 0;
end
params.customMerge = P.customMerge;


if ~isfield(P, 'performCorr')
    P.performCorr = 0;
end
params.performCorr = P.performCorr;

if params.performCorr
    if ~isfield(P, 'saveHighlight')
        P.saveHighlight = 0;
    end
    params.saveHighlight = P.saveHighlight;
    
    if ~isfield(P,'padCorrelation')
        P.padCorrelation = 1;
    end
    params.padCorrelation = P.padCorrelation;
end


if isfield(P, 'brightenNorming') && P.brightenNorming
    params.brightenNorming = 1;
    
    if ~isfield(P, 'brightConstant')
        P.brightConstant = 3000;
    end
    
    if ~isfield(P, 'brightGaussSDs')
        % Guassian width in pixels: 
%         50 * [512 402] ./ [710 690]; % 50 microns; scan voltage: 4.9 (MTK)
%         30 * [512 512] ./ [585 585]; % 30 microns; scan voltage: 4 (FN)
        P.brightGaussSDs = 50 * [512 402] ./ [710 690]; % MTK: scan voltage: 4.9
%         P.brightGaussSDs = 30 * [512 512] ./ [585 585]; % FN: scan voltage: 4
        % Default is 50 um scaled to MTK scan settings
        % (The "50" is the 50 um, the second term is to convert microns to
        % pixels, in height and width, respectively)
    end
    
    params.brightGaussSDs = P.brightGaussSDs;
    params.brightConstant = P.brightConstant;

    if ischar(P.brightConstant) && strcmp(P.brightConstant, 'Quantile')
        if ~isfield(P, 'brightConstantQuantile')
            P.brightConstantQuantile = 0.5;
        end
        params.brightConstantQuantile = P.brightConstantQuantile;
    end
    
else
    params.brightenNorming = 0;
end


if ~isfield(P, 'performCorrEP') % Efty's version of our correlationThreshold code.
    P.performCorrEP = 0;
end
params.performCorrEP = P.performCorrEP;


if ~isfield(P, 'cnmfStep1Done') % Whether CNMF steps before deconvolution are already done.
    P.cnmfStep1Done = false;
end
params.cnmfStep1Done = P.cnmfStep1Done;


if ~isfield(P, 'pnevFileName')
    P.pnevFileName = '';
end
params.pnevFileName = P.pnevFileName;
if params.cnmfStep1Done && strcmp(params.pnevFileName, '')
    error('You need to specify pnevFileName.')
end


if ~isfield(P, 'tau')
    P.tau = 4; % std of gaussian kernel (size of neuron)
end
params.tau = P.tau; 


if ~isfield(P, 'makeMCMrepMovie')
    P.makeMCMrepMovie = 0;
end
params.makeMCMrepMovie = P.makeMCMrepMovie;


if params.saveGoodMovieStats==1 || ~isfield(P, 'onlyAveImg') % FN: P.onlyAveImg=1; all projection images are saved except for aveImage, and we want to set it.
    P.onlyAveImg = 0;
end
params.onlyAveImg = P.onlyAveImg;

if ~isfield(P, 'frsExclude') % frames not to use when computing the mask (related to motion correction).
    P.frsExclude = [];
end
params.frsExclude = P.frsExclude;


%% Save

if isfield(P, 'saveParams')
    saveParams = P.saveParams;
else
    saveParams = true;
end

if saveParams && ~isempty(outName)
    if ~exist('improcparams', 'dir')
        warning('Folder improcparams does not exist in the current directory; so params not saved!')
    else
        save(fullfile('improcparams', outName), 'params');
    end
end


%% Copy to workspace

assignin('base', 'params', params);
