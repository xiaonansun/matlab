function [A,b,C,P] = update_spatial_components(Y,C,f,A_,P,options)

% update spatial footprints and background through Basis Pursuit Denoising
% for each pixel i solve the problem 
%   [A(i,:),b(i)] = argmin sum(A(i,:))
%       subject to || Y(i,:) - A(i,:)*C + b(i)*f || <= sn(i)*sqrt(T);
% for each pixel the search is limited to a few spatial components

% INPUTS:
% Y:    raw data
% C:    temporal components
% f:    temporal background
% A_:   current estimate of spatial footprints (used for determining search locations only) 
% P:    dataset parameters (used for noise values and interpolated entries)

% options    parameter struct (for noise values and other parameters)

% OUTPUTS:
% A:    new estimate of spatial footprints
% b:    new estimate of spatial background
% C:    temporal components (updated only when spatial components are completely removed)

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
memmaped = isobject(Y);
if memmaped
    sizY = size(Y,'Y');
    d = prod(sizY(1:end-1));
    T = sizY(end);
else
    [d,T] = size(Y);
end
if nargin < 6 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); options.d1 = d1; else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); options.d2 = d2; else d2 = options.d2; end          % # of columns
if ~isfield(options,'d3') || isempty(options.d3); d3 = input('What is the total number of z-planes? \n'); options.d3 = d3; else d3 = options.d3; end          % # of columns
if ~isfield(options,'show_sum'); show_sum = 0; else show_sum = options.show_sum; end            % do some plotting while calculating footprints
if ~isfield(options,'interp'); Y_interp = sparse(d,T); else Y_interp = options.interp; end      % identify missing data
if ~isfield(options,'spatial_parallel'); spatial_parallel = ~isempty(which('parpool')); else spatial_parallel = options.spatial_parallel; end % use parallel toolbox if present
if ~isfield(options,'search_method'); method = []; else method = options.search_method; end     % search method for determining footprint of spatial components
if ~isfield(options,'tsub') || isempty(options.tsub); tsub = 1; else tsub = options.tsub; end  % downsample temporally to estimate A and b

if nargin < 2 || (isempty(A_) && isempty(C))  % at least either spatial or temporal components should be provided
    error('Not enough input arguments')
else
    if ~isempty(C); K = size(C,1); elseif islogical(A_); K = size(A_,2); else K = size(A_2,2) - options.nb; end
end

if nargin < 5 || isempty(P); P = preprocess_data(Y,1); end  % etsimate noise values if not present


%%
% FN: We start parpool here (and specify its numWorkers and numThreads) because determine_search_location uses parfor.
if spatial_parallel
    % FN : create parallel pool of determined size
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        if options.poolsize
            % if no pools exist, create one with NumWorkers equal options.poolsize
            poolobj = parpool(options.poolsize);
        else
            poolobj = parpool; % go with matlab's default numWorkers, not recommended when running jobs on the shared cluster.
        end
    else
        fprintf('parpool already exists.\n')
    end
    poolsize = poolobj.NumWorkers;
    fprintf('Parallel pool size = %i\n', poolsize)
    % FN: increase the IdleTimeout of parpool.
    poolobj.IdleTimeout = Inf;
    
    % MTK: if requested, limit threads per worker
    % FN modified to account for numWorkers in parallel pool.
    fprintf('maxNumCompThreads before parfor = %.1d\n', maxNumCompThreads) % FN check for num threads.
    if options.limit_threads
        origThreadCount = maxNumCompThreads; % FN: this should be what you have set at the beginning of demo_script_modif. We go back to this number at the end of parfor.
        % FN: outside parfor we want limit_thrads to be 16, inside parfor we want it to be 16/8=2 if asked for 8 workers.
        newThreadNums = ceil(options.limit_threads / poolsize);
        maxNumCompThreads(newThreadNums);
    end
    fprintf('maxNumCompThreads, reset for parfor = %.1d\n', maxNumCompThreads) % FN check for num threads.     
end


%%
if nargin < 4 || isempty(A_); 
    IND = ones(d,size(C,1)); 
else
    if islogical(A_)     % check if search locations have been provided, otherwise estimate them
        IND = A_;
        if isempty(C)    
            INDav = double(IND)/diag(sum(double(IND)));          
            px = (sum(IND,2)>0);
            f = mean(Y(~px,:));
            b = max(Y*f',0)/norm(f)^2;
            C = max(INDav'*Y - (INDav'*b)*f,0);
        end
    else
        fprintf('Determine search locations...') % FN
        IND = determine_search_location(A_(:,1:K),method,options);
        fprintf('... Done.\n')
    end
end

Cf = [C;f];

if tsub ~= 1
    P.sn_ds = zeros(d,1);
    Ts = floor(T/tsub);
    Cf = squeeze(mean(reshape(Cf(:,1:Ts*tsub),[],tsub,Ts),2));        
    f = Cf(end-size(f,1)+1:end,:);
    T = Ts;
end 

options.sn = P.sn;
if ~memmaped
    Y(P.mis_entries) = NaN; % remove interpolated values
    if tsub~=1        
        Y = squeeze(mean(reshape(Y(:,1:Ts*tsub),[],tsub,Ts),2));    % downsample
        if ~isfield(P,'sn_ds');
            [P_ds,Y] = preprocess_data(Y);            
            P.sn_ds = P_ds.sn;        
        end
        options.sn = P.sn_ds;
    end
end

fprintf('Update A ...\n') % FN   
if spatial_parallel         % solve BPDN problem for each pixel
    Nthr = max(20*maxNumCompThreads,round(d*T/2^24));
    Nthr = min(Nthr,round(d/1e3));
    siz_row = [floor(d/Nthr)*ones(Nthr-mod(d,Nthr),1);(floor(d/Nthr)+1)*ones(mod(d,Nthr),1)];
    indeces = [0;cumsum(siz_row)];
    Yf = cell(Nthr,1);
    A = spalloc(d,size(Cf,1),nnz(IND)+size(Cf,1)*d);
    for nthr = 1:Nthr     
        if memmaped
            Ytemp = double(Y.Yr(indeces(nthr)+1:indeces(nthr+1),:));
            if tsub ~= 1;
                Ytemp = squeeze(mean(reshape(Ytemp(:,1:Ts*tsub),[],tsub,Ts),2));
                if ~isfield(P,'sn_ds')
                    [P_ds,Ytemp] = preprocess_data(Ytemp);
                    P.sn_ds(indeces(nthr)+1:indeces(nthr+1)) = P_ds.sn;
                end
                sn_temp = P.sn_ds(indeces(nthr)+1:indeces(nthr+1));
            else
                sn_temp = P.sn(indeces(nthr)+1:indeces(nthr+1));
            end
        else
            Ytemp = Y(indeces(nthr)+1:indeces(nthr+1),:);
            sn_temp = options.sn(indeces(nthr)+1:indeces(nthr+1));
        end
        IND_temp = IND(indeces(nthr)+1:indeces(nthr+1),:);
        Atemp = spalloc(siz_row(nthr),size(Cf,1),nnz(IND_temp));
        Yf{nthr} = Ytemp*f'; 
        
        parfor px = 1:siz_row(nthr)
            fn = ~isnan(Ytemp(px,:));       % identify missing data
            ind = find(IND_temp(px,:));
            if ~isempty(ind);
                ind2 = [ind,K+(1:size(f,1))];
                [~, ~, a, ~] = lars_regression_noise(Ytemp(px,fn)', Cf(ind2,fn)', 1, sn_temp(px)^2*T);
                a_sparse = sparse(1,ind2,double(a'));
                Atemp(px,:) = a_sparse';
            end
        end
        if mod(nthr,50) == 0
            fprintf('%2.1f%% of pixels completed \n', indeces(nthr+1)*100/d);
        end
        %Acell{nthr} = Atemp;
        A(indeces(nthr)+1:indeces(nthr+1),:) = Atemp;
    end
    
    % MTK: restore max number of threads if needed
    if options.limit_threads
        maxNumCompThreads(origThreadCount);
    end
    fprintf('maxNumCompThreads after parfor = %.1d\n', maxNumCompThreads) % FN check for num threads.
    
    %A = cell2mat(Acell);
    Yf = cell2mat(Yf);
else
    A = [zeros(d,K),zeros(d,size(f,1))];
    sA = zeros(d1,d2,d3);
    Yf = Y*f';
    for px = 1:d   % estimate spatial components
        fn = ~isnan(Y(px,:));       % identify missing data
        ind = find(IND(px,:));
        if ~isempty(ind);
            ind2 = [ind,K+(1:size(f,1))];
            [~, ~, a, ~] = lars_regression_noise(Y(px,fn)', Cf(ind2,fn)', 1, options.sn(px)^2*T);
            A(px,ind2) = a';
            sA(px) = sum(a);
        end
        if show_sum
            if mod(px,d1) == 0;
               figure(20); imagesc(max(sA,[],3)); axis square;  
               title(sprintf('Sum of spatial components (%i out of %i columns done)',round(px/d1),d2)); drawnow;
            end
        end
    end
end

A(isnan(A))=0;
A = sparse(A);
A = threshold_components(A,options);  % post-processing of components

fprintf('Updated spatial components \n');

ff = find(sum(A(:,1:K))==0);           % remove empty components
if ~isempty(ff)
    K = K - length(ff);
    A(:,ff) = [];
    C(ff,:) = [];
    Cf(ff,:) = [];
end

b = double(max((Yf - A(:,1:K)*double(Cf(1:K,:)*f'))/(f*f'),0));
A = A(:,1:K);