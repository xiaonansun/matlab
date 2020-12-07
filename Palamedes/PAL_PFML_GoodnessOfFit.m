%
%PAL_PFML_GoodnessOfFit     Determine Goodness-of-Fit of psychometric
%   function (PF) using method in Wichmann & Hill, 2001, Perception &
%   Psychophysics, 63, 1293-1313.
%
%Syntax:    [Dev pDev DevSim converged] = ...
%           PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
%           paramsValues, paramsFree, B, PF,{optional arguments});
%
%Input:
%   'StimLevels': vector containing stimulus levels used.
%
%   'NumPos': vector containing for each of the entries of 'StimLevels' the 
%       number of trials a positive response (e.g., 'yes' or 'correct') was
%       given.
%
%   'OutOfNum': vector containing for each of the entries of 'StimLevels' 
%       the total number of trials.
%
%   'paramsValues': 1x4 vector containing parametervalues [threshold slope 
%       guess-rate lapse-rate]. PAL_PFML_Fit might be used to obtain
%       best-fitting parameter values.
%
%   'paramsFree': 1x4 vector coding which of the four parameters of the PF 
%       (in the order: [threshold slope guess-rate lapse-rate]) are free 
%       parameters and which are fixed parameters (1: free, 0: fixed, ).
%
%   'B': number of bootstrap simulations to perform.
%
%   'PF': psychometric function used in fit. Passed as an inline function.
%       Options include:    
%           @PAL_Logistic
%           @PAL_Weibull
%           @PAL_Gumbel (i.e., log-Weibull)
%           @PAL_Quick
%           @PAL_logQuick
%           @PAL_CumulativeNormal
%           @PAL_Gumbel
%           @PAL_HyperbolicSecant
%
%Output: 
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       psychometric function to fit of saturated model)
%
%   'pDev': proportion of the B Deviance values from simulations that were
%       greater than Deviance value of data. The greater the value of pDev,
%       the better the fit.
%
%   'DevSim': vector containing all B simulated Deviance values.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   PAL_PFML_GoodnessOfFit will generate a warning if not all 
%       simulations were fit succesfully.
%
%   PAL_PFML_GoodnessOfFit will accept a few optional arguments:
%
%       Use 'searchGrid' argument to define a 4D parameter 
%       grid through which to perform a brute-force search for initial 
%       guesses (performed by PAL_PFML_BruteForceFit) to be used during 
%       fitting procedure. Structure should have fields .alpha, .beta, 
%       .gamma, and .lambda. Each should list parameter values to be 
%       included in brute force search. Fields for fixed parameters should 
%       be scalars equal to the fixed parameter value. Note that all fields 
%       may be scalars in which case no brute-force search will precede the 
%       iterative parameter search. For more information, see example 
%       below. Note that choices made here may have a large effect on 
%       processing time and memory usage. Note that usage of the 
%       'searchGrid' option will generally reduce processing time ((because
%       the serial iterative Nelder-Mead search will be shorter). Usage of 
%       the 'searchGrid' option will also reduce the chances of a fit 
%       failing (see: www.palamedestoolbox.org/understandingfitting.html).
%       In some future version of Palamedes the 'searchGrid' argument will 
%       be required. 
%
%       If highest entry in 'StimLevels' is so high that it can be 
%       assumed that errors observed there can be due only to lapses, use 
%       'lapseFit' argument to specify alternative fitting scheme. Options: 
%       'nAPLE' (default), 'iAPLE', and 'jAPLE'. Type help 
%       PAl_PFML_FitMultiple for more information.
%
%       The guess rate and lapse rate parameter can be constrained to be 
%       equal, as would be appropriate, for example, in a bistable percept 
%       task. To accomplish this, use optional argument 'gammaEQlambda', 
%       followed by a 1. Both the guess rate and lapse rate parameters will 
%       be fit according to options set for the lapse rate parameter. Entry
%       for guess rate in 'searchGrid' needs to be made but will be 
%       ignored.
%
%       Options 'maxTries' and 'rangeTries' were operational up to and
%       including Palamedes version 1.7.0. Usage of these arguments will
%       result in a warning stating that the arguments will be ignored and
%       will encourage users to utilize the 'searchGrid' options (see
%       above), which is a much superior manner in which to avoid failed
%       fits.
%
%       User may constrain the lapse rate to fall within limited range 
%       using the optional argument 'lapseLimits', followed by a two-
%       element vector containing lower and upper limit respectively. See 
%       full example below.
%
%       User may constrain the guess rate to fall within limited range 
%       using the optional argument 'guessLimits', followed by a two-
%       element vector containing lower and upper limit respectively.
%
%   PAL_PFML_BootstrapNonParametric uses Nelder-Mead Simplex method to find 
%   the maximum in the likelihood function. The default search options may 
%   be changed by using the optional argument 'SearchOptions' followed by 
%   an options structure created using options = PAL_minimize('options'), 
%   then modified. See example of usage in PAL_PFML_Fit. For more 
%   information type PAL_minimize('options','help').
%
%       If the likelihood function contains a global maximum (and assuming 
%       that an appropriate search grid is used), the fitting procedure 
%       will find the global maximum. However, sometimes the likelihood 
%       function does not contain a global maximum. This situation would 
%       occasionally result in 'false' convergences (procedure claims to 
%       have found global maximum but has not). Starting in Palamedes 
%       version 1.10.0 fitting procedures will, by default, check whether a 
%       fit corresponds to a global maximum. If not, the procedure will 
%       identify the best-fitting function that can be asymptotically 
%       approached by the PF (this will be either a step function or a 
%       constant function). Even though this poses problems for assigning 
%       finite values to maximum-likelihood parameter estimates, the log-
%       likelihood value for such a function is finite and can be 
%       asymptotically approached by a psychometric function. Thus,
%       the absence of a true maximum is not really a problem for 
%       determining Goodness-of-Fit based on simulated log-likelihood
%       values (see www.palamedestoolbox.understandingfitting.html for
%       more information). Even though there is really no reason to do so, 
%       user may override the default behavior of checking whether true 
%       maximum exists by using the 'checkLimits' option followed by 0 (or 
%       logical false). It is important to note that the lack of a true
%       maximum in the likelihood may be the result of 'Too much model, too 
%       little data' (i.e., one is trying to estimate parameters that the 
%       data do not contain enough information on). A better solution
%       is to either get more data or fit a simpler model. Another option
%       is to use a Bayesian fitting approach.
%
%Example:
%
%   PF = @PAL_Logistic;
%   StimLevels = [-3:1:3];
%   NumPos = [55 55 66 75 91 94 97];    %observer data
%   OutOfNum = 100.*ones(size(StimLevels));
%   searchGrid.alpha = [-1:.1:1];    %structure defining grid to
%   searchGrid.beta = 10.^[-1:.1:2]; %search for initial values
%   searchGrid.gamma = .5;
%   searchGrid.lambda = [0:.005:.03];
%   paramsFree = [1 1 0 1];
%
%   %First fit data:
%
%   paramsValues = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, ...
%       searchGrid, paramsFree, PF,'lapseLimits',[0 .03]);
%
%   %Determine Goodness-Of-Fit
%
%   B = 400;
%
%   [Dev pDev DevSim converged] = ...
%       PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid, ...
%       'lapseLimits',[0 .03]);
%
% Introduced: Palemedes version 1.0.0 (NP)
% Modified: Palamedes version 1.0.2, 1.1.0, 1.2.0, 1.3.9, 1.3.1, 1.4.0,
%   1.6.3, 1.8.1, 1.9.1 (see History.m)

function [Dev, pDev, DevSim, converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, paramsValues, paramsFree, B, PF, varargin)

searchGrid = paramsValues;

options = [];
lapseLimits = [];
guessLimits = [];
lapseFit = 'default';
gammaEQlambda = logical(false);
checkLimits = paramsFree(2);

converged = zeros(B,1);
DevSim = zeros(B,1);

mTrTflag = false;

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strcmpi(varargin{n}, 'SearchOptions')
            options = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'maxTries',4)
            mTrTflag = true;
            valid = 1;
        end
        if strncmpi(varargin{n}, 'rangeTries',6)
            mTrTflag = true;
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseLimits',6)
            if paramsFree(4) == 1
                lapseLimits = varargin{n+1};
            else
                warning('PALAMEDES:invalidOption','Lapse rate is not a free parameter: ''LapseLimits'' argument ignored');
            end
            valid = 1;
        end
        if strncmpi(varargin{n}, 'guessLimits',6)
            if paramsFree(3) == 1
                guessLimits = varargin{n+1};
            else
                warning('PALAMEDES:invalidOption','Guess rate is not a free parameter: ''GuessLimits'' argument ignored');
            end
            valid = 1;
        end
        if strncmpi(varargin{n}, 'searchGrid',8)
            searchGrid = varargin{n+1};
            valid = 1;
        end        
        if strncmpi(varargin{n}, 'lapseFit',6)
            if paramsFree(4) == 0  && ~strncmp(varargin{n+1},'def',3)
                warning('PALAMEDES:invalidOption','Lapse rate is not a free parameter: ''LapseFit'' argument ignored');   
            else
                if strncmpi(varargin{n+1}, 'nAPLE',5) || strncmpi(varargin{n+1}, 'jAPLE',5) || strncmpi(varargin{n+1}, 'default',5)
                    lapseFit = varargin{n+1};
                else 
                    if strncmpi(varargin{n+1}, 'iAPLE',5)
                        warning('PALAMEDES:invalidOption','iAPLE fitting no longer supported, using jAPLE fitting instead (iAPLE instead of jAPLE fitting is hard to justify anyway).');
                    else
                        warning('PALAMEDES:invalidOption','%s is not a valid option for ''lapseFit''. ignored', varargin{n+1});   
                    end
                end                
            end 
            valid = 1;
        end
        if strncmpi(varargin{n}, 'gammaEQlambda',6)
            gammaEQlambda = logical(varargin{n+1});
            if gammaEQlambda                
                if paramsValues(3) ~= paramsValues(4)
                    paramsValues(3) = paramsValues(4);
                    warning('PALAMEDES:invalidOption','Generating gamma value changed to %s in order to match lapse value.',num2str(paramsValues(3)));
                end                
            valid = 1;
            end
        end  
        if strncmpi(varargin{n}, 'checkLimits',6)
            checkLimits = varargin{n+1};
            valid = 1;
        end        
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end
    end            
end

if mTrTflag
    message = ['''maxTries'' and ''rangeTries'' options are obsolete and '];
    message = [message 'will be ignored. Use ''searchGrid'' option instead. See: '];
    message = [message 'www.palamedestoolbox.org/understandingfitting.html '];
    message = [message 'for more information.'];
    warning('PALAMEDES:invalidOption',message);
end

if ~isempty(guessLimits) && gammaEQlambda
    warning('PALAMEDES:invalidOption','Guess rate is constrained to equal lapse rate: ''guessLimits'' argument ignored');
    guessLimits = [];
end

if ~isstruct(searchGrid)
    message = ['Option to use generating parameter values as initial '];    
    message = [message 'guesses in fitting procedure will be removed in '];
    message = [message 'some future version of Palamedes. Instead use '];
    message = [message 'optional argument ''searchGrid'' to pass a '];
    message = [message 'structure defining a 4D parameter space through '];
    message = [message 'which to search for initial search values using a '];
    message = [message 'brute force search. Type help '];
    message = [message 'PAL_PFML_GoodnessOfFit for more information.'];
    warning('PALAMEDES:useSearchGrid',message);
end

[StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum);

negLLCon = PAL_PFML_negLL([], paramsValues, [0 0 0 0], StimLevels, NumPos, OutOfNum, PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
negLLAug = PAL_PFML_negLLNonParametric(NumPos, OutOfNum);

Dev = 2*(negLLCon-negLLAug);

if isstruct(searchGrid)
    
    if gammaEQlambda
        searchGrid.gamma = 0;
    end
    
    [paramsGrid.alpha, paramsGrid.beta, paramsGrid.gamma, paramsGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda);

    singletonDim = uint16(size(paramsGrid.alpha) == 1);    
    
    [paramsGrid.alpha] = squeeze(paramsGrid.alpha);
    [paramsGrid.beta] = squeeze(paramsGrid.beta);
    [paramsGrid.gamma] = squeeze(paramsGrid.gamma);
    [paramsGrid.lambda] = squeeze(paramsGrid.lambda);

    if gammaEQlambda
        paramsGrid.gamma = paramsGrid.lambda;
    end
    
    for level = 1:length(StimLevels)
        logpcorrect(level,:,:,:,:) = log(PF(paramsGrid,StimLevels(1,level)));
        logpincorrect(level,:,:,:,:) = log(1-PF(paramsGrid,StimLevels(1,level)));
    end
else
    paramsGuess = paramsValues;
    if gammaEQlambda
        paramsGuess(3) = paramsGuess(4);
        paramsFree(3) = 0;
    end        
end

for b = 1:B
    NumPosSim = PAL_PF_SimulateObserverParametric(paramsValues, StimLevels, OutOfNum, PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    
    if isstruct(searchGrid)
        LLspace = zeros(size(paramsGrid.alpha,1),size(paramsGrid.alpha,2),size(paramsGrid.alpha,3),size(paramsGrid.alpha,4));

        if size(LLspace,2) == 1 && ndims(LLspace) == 2
            LLspace = LLspace';
        end

        switch lower(lapseFit(1:3))
            case {'nap','def'}
                for level = 1:length(StimLevels)
                   LLspace = LLspace + NumPosSim(level).*squeeze(logpcorrect(level,:,:,:,:))+(OutOfNum(level)-NumPosSim(level)).*squeeze(logpincorrect(level,:,:,:,:));
                end
            case {'jap','iap'}
                len = length(NumPosSim);                            
                LLspace = log((1-paramsGrid.lambda).^NumPosSim(len))+log(paramsGrid.lambda.^(OutOfNum(len)-NumPosSim(len))); %0*log(0) evaluates to NaN, log(0.^0) does not
                if gammaEQlambda
                    LLspace = LLspace + log(paramsGrid.lambda.^NumPosSim(1)) + log((1-paramsGrid.lambda).^(OutOfNum(1)-NumPosSim(1)));
                end
                for level = 1+gammaEQlambda:len-1
                    LLspace = LLspace + NumPosSim(level).*squeeze(logpcorrect(level,:,:,:,:))+(OutOfNum(level)-NumPosSim(level)).*squeeze(logpincorrect(level,:,:,:,:));
                end    
        end
        
        if isvector(LLspace)
            [maxim, Itemp] = max(LLspace);
        else
            if strncmpi(lapseFit,'iap',3)
                [trash, lapseIndex] = min(abs(searchGrid.lambda-(1-NumPosSim(len)/OutOfNum(len))));
                LLspace = shiftdim(LLspace,length(size(LLspace))-1);
                [maxim, Itemp] = PAL_findMax(LLspace(lapseIndex,:,:,:));
                Itemp = circshift(Itemp',length(size(LLspace))-1)';
            else
                [maxim, Itemp] = PAL_findMax(LLspace);
            end
        end
        
        I = ones(1,4);
         
        I(singletonDim == 0) = Itemp;                
            
        paramsGuess = [searchGrid.alpha(I(1)) searchGrid.beta(I(2)) searchGrid.gamma(I(3)) searchGrid.lambda(I(4))];
        if strncmpi(lapseFit,'iap',3)
            paramsGuess(4) = 1-NumPosSim(len)/OutOfNum(len);
        end

    end    
    
    [paramsValuesSim, trash, converged(b)] = PAL_PFML_Fit(StimLevels, NumPosSim, OutOfNum, paramsGuess, paramsFree, PF, 'SearchOptions', options,'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda,'checkLimits',checkLimits);

    negLLConSim = PAL_PFML_negLL([], paramsValuesSim, [0 0 0 0], StimLevels, NumPosSim, OutOfNum, PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    negLLAugSim = PAL_PFML_negLLNonParametric(NumPosSim, OutOfNum);
    DevSim(b) = 2*(negLLConSim-negLLAugSim);
    if converged(b) ~= 1
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end

end

if any(converged ~= 1)
    exitflag = 0;
    message = [char(10), char(10), 'Only %s of %s simulations converged successfully. ', char(10)];
    message = [message, 'Only fits for which entry in ''converged'' equals 1 were successful (i.e.,', char(10)];
    message = [message, 'a global maximum in likelihood function was found).', char(10)];
    message = [message, '%s simulations resulted in scenario -1', char(10)];
    message = [message, '%s simulations resulted in scenario -2', char(10)];
    message = [message, '%s simulations resulted in scenario -3', char(10)];
    message = [message, 'However, this is not really a problem for interpretation of Dev or pDev', char(10)];
    message = [message, 'values (all simulated log-likelihood values can be approached asymptotically ', char(10)];
    message = [message, 'by the model).', char(10)];
    message = [message, 'Nevertheless, particularly elegant it is not.', char(10)];
    message = [message, 'See www.palamedestoolbox.org/understandingfitting.html for more '];
    message = [message, 'information', char(10)];
    message = [message, 'Some possible solutions:', char(10)];
    message = [message, '-Adjust the searchGrid option.', char(10)];
    message = [message, '-Fix one or more parameters (lapse rate is a good first choice)', char(10)];
    message = [message, '-Collect more data', char(10)];
    if ~checkLimits
        message = [message, 'Use ''checkLimits'' option in this function to get a better', char(10)];    
        message = [message, 'idea of why fits failed (type ''help PAL_PFML_GoodnessOfFit'').', char(10)];    
    end
    warning('PALAMEDES:convergeFail',message,int2str(sum(converged == 1)), int2str(B),int2str(sum(converged == -1)),int2str(sum(converged == -2)),int2str(sum(converged == -3)));
else
    exitflag = 1;    
end

pDev = length(DevSim(DevSim>Dev))/B;