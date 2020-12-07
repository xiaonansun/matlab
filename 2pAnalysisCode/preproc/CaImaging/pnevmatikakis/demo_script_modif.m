function [A, C, S, C_df, Df, b, f, srt, srt_val, Ain, options, P, merging_vars, YrA, Yr, Cin, bin, fin, nA,AY,AA] = demo_script_modif(Y, pnev_inputParams)
% [A, C, S, C_df, S_df, Df, b, f, srt, options, P] = demo_script_modif(Y, pnev_inputParams);
%
% INPUTS:
%     Y                 % movie (height x width x frames)
%     pnev_inputParams  % structure of input parameters with the following fields:
%         K                      % number of components to be found
%         temp_sub               % temporal subsampling for greedy initiation, set to 1 for no down sampling.
%         space_sub              % spatial subsampling for greedy initiation, set to 1 for no down sampling.
%         tau                    % std of gaussian kernel (size of neuron)
%         p                      % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
%         merge_thr              % merging threshold
%         deconv_method          % activity deconvolution method
%         temp_iter              % number of block-coordinate descent steps
%         fudge_factor           % bias correction for AR coefficients
%         finalRoundMCMC         % do a final round of MCMC method (if false, after merging 2 iterations of const foopsi will be done. If true, after merging 1 iter of const foopsi and 1 iter of MCMC will be done.)
%         doPlots                % make some figures and a movie.
%         parallelTempUpdate     % do parallel temporal updating.
%         save4debug             % save Eftychios's variables (eg A,C,etc) after each step for debug purposes.
%         search_dist            % search distance when updating spatial components.
%
%
% OUTPUTS:
% Note: outputs A,C,S are normalized to nA (ie norm of A).
%     A;          % spatial components of neurons
%     C;          % temporal components of neurons
%     S;          % spike counts
%     C_df;       % temporal components of neurons and background normalized by Df
%     S_df;       % spike counts of neurons normalized by Df
%     Df;         % background for each component to normalize the filtered raw data
%     b;          % spatial components of backgrounds
%     f;          % temporal components of backgrounds
%     srt         % index of ROIs before ordering
%     Ain         % spatial components of neurons after initialization
%     options;    % options for model fitting
%     P;          % some estimated parameters
%     merging_vars % some vars related to merging; allows offline assessment if desired.
%
% based on Eftychios's demo_script (repository ca_source_extraction V0.2.1)
% use processCaImagingMCPnev.m to set INPUT.
% FN (Jan 11 2016)



%% Set maxNumCompThreads

% FN: Choose limit_threads similar to -pe threads, ie ~16.
fprintf('Original maxNumCompThreads = %.1d\n', maxNumCompThreads)
if pnev_inputParams.limit_threads
    origThreadCount0 = maxNumCompThreads;
    maxNumCompThreads(pnev_inputParams.limit_threads);
    fprintf('maxNumCompThreads after initial reset = %.1d\n', maxNumCompThreads)
end


%% Subtract out the min value from the movie and convert it to single

t1 = tic;
if isempty(pnev_inputParams.maxFrsForMinPsn)
    Y = Y - min(Y(:));
else
    Y = Y - min(reshape(Y(:,:, 1:min(pnev_inputParams.maxFrsForMinPsn, size(Y,3))), [], 1));
end

if ~isa(Y,'single');    Y = single(Y);  end         % convert to single % FN: As of 7/29/16. Before that it was double.

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

fprintf('\nConverting to single and subtracting min took %0.1f s\n\n', toc(t1));


%% Set options

K = pnev_inputParams.K;             % number of components to be found
tau = pnev_inputParams.tau;         % std of gaussian kernel (size of neuron); Eftychios's default=4;
p = pnev_inputParams.p;             % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay); Eftychios's default=2;
merge_thr = pnev_inputParams.merge_thr; % Eftychios's default=0.8;
deconv_meth = pnev_inputParams.deconv_method;
temp_iter = pnev_inputParams.temp_iter;
fudge_fact = pnev_inputParams.fudge_factor;
temp_sub = pnev_inputParams.temp_sub;
space_sub = pnev_inputParams.space_sub;
% finalRoundMCMC = pnev_inputParams.finalRoundMCMC;
doPlots = pnev_inputParams.doPlots;
% parallelTempUpdate = pnev_inputParams.parallelTempUpdate; % FN commented.
% In V0.3.3 options.temporal_parallel takes care of it and it will be by
% default true if parallel cores exist on the machine.
temp_par = pnev_inputParams.temporal_parallel;
save4debug = pnev_inputParams.save4debug;
search_dist = pnev_inputParams.search_dist;
init_method = pnev_inputParams.init_method;
orderROI_extractDf = pnev_inputParams.orderROI_extractDf;
noise_norm = pnev_inputParams.noise_norm;
noise_norm_prctile = pnev_inputParams.noise_norm_prctile;
bas_nonneg = pnev_inputParams.bas_nonneg;

options = CNMFSetParms(...
    'd1',d1,'d2',d2,...                                   % dimensions of datasets
    'search_method','ellipse','dist',search_dist,...      % search locations when updating spatial components
    'deconv_method',deconv_meth,...                       % activity deconvolution method
    'temporal_iter',temp_iter,...                         % number of block-coordinate descent steps
    'fudge_factor',fudge_fact,...                         % bias correction for AR coefficients
    'merge_thr',merge_thr,...                             % merging threshold
    'gSig',tau,...
    'tsub', temp_sub,...                                  % temporal subsampling for greedy initiation
    'ssub', space_sub, ...                                % spatial subsampling for greedy initiation
    'init_method', init_method, ...                       % initialization method
    'temporal_parallel', temp_par, ...                    % flag for parallel updating of temporal components (default: true if present)
    'noise_norm', noise_norm, ...                         % whether to do normalization by noise estimate prior to initialization. (default: true)
    'noise_norm_prctile', noise_norm_prctile, ...         % minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
    'bas_nonneg', bas_nonneg ...                          % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
    );

% limit_threads is an option added by MTK to limit the number of cores used
% per parallel thread. Matlab, when running many workers on a single
% machine, uses MaxNumCompThreads for each worker (it seems). This is a
% terrible idea, and leads to poor performance as num_workers * num_cores
% threads are swapped in and out of num_cores cores.
%
% limit_threads lets you set the number of cores to use per worker. Likely
% values are 2-4. 0 disables limiting (default)
if isfield(pnev_inputParams, 'limit_threads')
    options.limit_threads = pnev_inputParams.limit_threads;
else
    options.limit_threads = 0;
end

% FN added the following optional parameter to allow the control of number
% of parallel workers when using parfor loops in update_spatial and
% temporal components.
if isfield(pnev_inputParams, 'poolsize')
    options.poolsize = pnev_inputParams.poolsize;
else
    options.poolsize = 0;
end

if isfield(pnev_inputParams, 'normingMedianImage')
    options.normingMedianImage = pnev_inputParams.normingMedianImage;
end

% FN: Below is commented because in the version released at the end of June
% 2016, Eftychios has added this feature to get_noise_fft using
% options.max_timesteps.
% FN added the option below to help with speed: when Y has many frames, we
% only use the 1st maxFrsForMinPsn frames to compute min(Y(:)) and estimate
% noies (P.sn).
% options.maxFrsForMinPsn = pnev_inputParams.maxFrsForMinPsn;


%% Data pre-processing

t1 = tic;

[P,Y] = preprocess_data(Y,p,options);

fprintf('\npreprocess_data took %0.1f s\n\n', toc(t1));


%% fast initialization of spatial components using greedyROI and HALS

% init_test
% options.normingMedianImage = pnev_inputParams.normingMedianImage;

t1 = tic;

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize % FN: providing P as input means that Y will be normalized to P.sn before initialization unless options.noise_norm is set to 0.

fprintf('\ninitialize_components took %0.1f s\n\n', toc(t1));


if doPlots
    % display centers of found components
    Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
    figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;
end


%% manually refine components (optional)

refine_components = false;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end


%% update spatial components

fprintf('Updating spatial components started.\n')

t1 = tic;

Yr = reshape(Y,d,T);
clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

fprintf('\nFirst update_spatial_components took %0.1f s\n\n', toc(t1));

if save4debug
    nowStr = datestr(now, 'yymmdd-HHMMSS');
    save(['Eft_aftInitAndSpUp-', nowStr], 'A','b','Cin','fin','P','options')
end


%% update temporal components

fprintf('Updating temporal components started.\n')

t1 = tic;

if pnev_inputParams.doMerging
    P.p = 0;    % set AR temporarily to zero for speed. The original value will be restored after merging.
end
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

fprintf('\nFirst update_temporal_components took %0.1f s\n\n', toc(t1));


t1 = tic;
complexTol = 1e-10;
[A, C, S, P, YrA] = removeComplexUnits(A, C, S, f, complexTol, P, YrA);

fprintf('\nRemoving complex units took %0.1f s\n\n', toc(t1));


if save4debug
    save(['Eft_preMerge-', nowStr], 'A', 'b', 'C', 'f', 'S', 'P', 'options')
end


%% merge found components

if pnev_inputParams.doMerging

    fprintf('Merging components started.\n')
    t1 = tic;
    
    if pnev_inputParams.customMerge 
        % use Farzaneh's function mergeROIs_set to define merged_ROIs.
%         th_dist_overlap_corr = [4.8 .5 .8]; % thresholds for COM distance, average mask overlap, and temporal correlation to find merged_ROIs.
        th_dist_overlap_Ccorr_Acorr = [4.8 .5 .75 .2]; %[4.8 .5 .8 .2];
        onlyCOM = 1; %0; % if 0, both COM dist and fraction overlap will be used.
        doEvaluate = 0;
        merged_ROIs = mergeROIs_set([], A, C, d1, d2, th_dist_overlap_Ccorr_Acorr, onlyCOM, doEvaluate, 0);        
    %     [Am,Cm,K_m,~,P,Sm] = merge_components_again(merged_ROIs,A,b,C,f,P,S,options);
        [Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options,merged_ROIs);
    else
        
        [Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);    
    end
    
    fprintf('\nMerging components took %0.1f s\n\n', toc(t1));

    
    display_merging = 1; % flag for displaying merging example
    if doPlots && display_merging && ~isempty(merged_ROIs)
        for i = 1:length(merged_ROIs); % randi(length(merged_ROIs));
            ln = length(merged_ROIs{i});
            figure;
            set(gcf,'Position',[300,300,(ln+2)*300,300]);
            for j = 1:ln
                subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2));
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
            end
            subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
            title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight;
            subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))');
            hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
            drawnow;
        end
    end
    
    if pnev_inputParams.save_merging_vars
        merging_vars.merged_ROIs = merged_ROIs;
        %     merging_vars.Am = Am;
        %     merging_vars.Cm = Cm;
        %     merging_vars.K_m = K_m;
%         merging_vars.A = A;
%         merging_vars.C = C;
    else
        merging_vars = [];
    end
    
    
    %% repeat: update_spatial_components
    
    P.p = p;    % restore AR value
    
    fprintf('Repeated updating of spatial components started.\n')
    t1 = tic;
    [A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,P,options);
    
    fprintf('\nSecond update_spatial_components took %0.1f s\n\n', toc(t1));
    
    if save4debug
        save(['Eft_aftMergeAndSpUp-', nowStr], 'A2', 'b2', 'Cm', 'f', 'Sm', 'P', 'options')
    end
    
    
    %% repeat: update_temporal_components
    
    fprintf('Repeated updating of temporal components started.\n')
    
    if ~P.p || ~strcmp(deconv_meth, 'MCMC') % ~finalRoundMCMC
        
        t1 = tic;
        [C2,f2,P,S2,YrA,AY,AA] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);
        fprintf('\nSecond update_temporal_components took %0.1f s\n\n', toc(t1));
        
        t1 = tic;
        [A2, C2, S2, P, YrA,AY,AA] = removeComplexUnits(A2, C2, S2, f2, complexTol, P, YrA,AY,AA);
        fprintf('\nRemoving complex units took %0.1f s\n\n', toc(t1));
        
        % save Cm and C2 for debugging
%         nowStr = datestr(now, 'yymmdd-HHMMSS');
%         save(['Eft_aftMergeAndSpTpUp-', nowStr], 'C2', 'Cm')
        
    else % FN: do 1 round const foopsi and then another round MCMC
        
        options.temporal_iter = 1;
        t1 = tic;
        [C2,f2,P,S2,YrA] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);
        fprintf('\nSecond update_temporal_components (1 round) took %0.1f s\n\n', toc(t1));
        
        [A2, C2, S2, P, YrA] = removeComplexUnits(A2, C2, S2, f2, complexTol, P, YrA);
        
        if save4debug
            save(['Eft_preMCMC-', nowStr], 'A2', 'b2', 'C2', 'f2', 'S2', 'P', 'options')
        end
        
        
%         options.deconv_method = 'MCMC';
        options.MCMC_B = pnev_inputParams.MCMC_B;
        options.MCMC_Nsamples = pnev_inputParams.MCMC_Nsamples;
        options.MCMC_prec = pnev_inputParams.MCMC_prec;
        fprintf('Final MCMC updating of temporal components started.\n')

        t1 = tic;
        [C2,f2,P,S2,YrA] = update_temporal_components(Yr,A2,b2,C2,f2,P,options);
        fprintf('\nFinal MCMC updating of temporal components took %0.1f s\n\n', toc(t1));
        
        [A2, C2, S2, P, YrA] = removeComplexUnits(A2, C2, S2, f2, complexTol, P, YrA);
        
    end
    
else
    A2 = A;
    C2 = C;
    S2 = S;
    b2 = b;
    f2 = f;
    merging_vars = [];
end


%%
if orderROI_extractDf
    
    % order ROIs. Also A, C S and P will be normalized by nA (ie norm of A).
    
    fprintf('Ordering ROIs...\n')
    t1 = tic;
%     [A_or, C_or, S_or, P, srt, srt_val, nA] = order_ROIs(A2,C2,S2,P, options);    % order components
    [A_or, C_or, S_or, YrA_or, P, srt, srt_val, nA] = order_ROIs(A2,C2,S2,YrA,P,options);
    fprintf('\norder_ROIs took %0.1f s\n\n', toc(t1));
%     YrA = YrA(srt,:);
    
    
    %% extract DF/F
    
    fprintf('Extracting DF/F...\n')
    t1 = tic;
    K_m = size(C_or, 1);
    [C_df, Df] = extract_DF_F(Yr, [A_or,b2], [C_or;f2], K_m+1); % extract DF/F values (optional) % FN moved it here so C_df and S_df are also ordered.
%     [C_df,Df,S_df] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],S_or,K_m+1);
    fprintf('\nextract_DF_F took %0.1f s\n\n', toc(t1));
    
  
else
    A_or = A2;
    C_or = C2;
    S_or = S2;
    YrA_or = YrA;
    C_df = []; Df = []; srt = []; srt_val = []; nA = []; 

end


%% do some plotting

if doPlots
    contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
    figure;
    [Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),contour_threshold,1); % contour plot of spatial footprints
%     pause;
    %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
    %     view_components(Yr,A_or,C_or,b2,f2,Cn,options);         % display all components
    
    % display components
    plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options)
end


%% make movie

if doPlots
    make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)
end


%%

% Remember A_or, C_or and S_or are normalized to nA (ie norm of A); hence,
% A,C,S (outputs of demo_script_modif) will be normalized to nA (ie norm of
% A) as well.
A = A_or;
C = C_or;
S = S_or;
YrA = YrA_or;

b = b2;
f = f2;

if ~exist('YrA', 'var') % || ~pnev_inputParams.multiTrs
    YrA = [];
end

if ~exist('AY', 'var')
    AY = [];
end

if ~exist('AA', 'var')
    AA = [];
end

if ~pnev_inputParams.multiTrs
    Yr = [];
end


%% Reset maxNumCompThreads to its original value

if options.limit_threads
    maxNumCompThreads(origThreadCount0);
end
fprintf('maxNumCompThreads, final reset = %.1d\n', maxNumCompThreads) % FN check for num threads.

