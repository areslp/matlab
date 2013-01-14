%% This function is modified from Matlab Package SpaRSA

function [x,e,iter]= SolveSpaRSA_CBM(A,y,varargin)

% SpaRSA version 2.0, December 31, 2007
%
% This function solves the convex problem
%
% arg min_x = 0.5*|| y - A x ||_2^2 + lambda phi(x)
%
% using the SpaRSA algorithm, which is described in "Sparse Reconstruction
% by Separable Approximation" by S. Wright, R. Nowak, M. Figueiredo,
% IEEE Transactions on Signal Processing, 2009 (to appear).
%
% The algorithm is related GPSR (Figueiredo, Nowak, Wright) but does not
% rely on the conversion to QP form of the l1 norm, because it is not
% limited to being used with l1 regularization. Instead it forms a
% separable
% approximation to the first term of the objective, which has the form
%
%  d'*A'*(A x - y) + 0.5*alpha*d'*d
%
% where alpha is obtained from a BB formula. In a monotone variant, alpha is
% increased until we see a decreasein the original objective function over
% this step.
%
% -----------------------------------------------------------------------
% Copyright (2007): Mario Figueiredo, Robert Nowak, Stephen Wright
%
% GPSR is distributed under the terms
% of the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%
% Please check for the latest version of the code and paper at
% www.lx.it.pt/~mtf/SpaRSA
%
%  ===== Required inputs =============
%
%  y: 1D vector or 2D array (image) of observations
%
%  A: if y and x are both 1D vectors, A can be a
%     k*n (where k is the size of y and n the size of x)
%     matrix or a handle to a function that computes
%     products of the form A*v, for some vector v.
%     In any other case (if y and/or x are 2D arrays),
%     A has to be passed as a handle to a function which computes
%     products of the form A*x; another handle to a function
%     AT which computes products of the form A'*x is also required
%     in this case. The size of x is determined as the size
%     of the result of applying AT.
%
%  lambda: regularization parameter (scalar)
%
%  ===== Optional inputs =============
%
%
%  'AT'    = function handle for the function that implements
%            the multiplication by the conjugate of A, when A
%            is a function handle. If A is an array, AT is ignored.
%
%  'Psi'   = handle to the denoising function, that is, to a function
%            that computes the solution of the densoing probelm
%            corresponding to the desired regularizer. That is,
%            Psi(y,lambda) = arg min_x (1/2)*(x - y)^2 + lambda phi(x).
%            Default: in the absence of any Phi given by the user,
%            it is assumed that phi(x) = ||x||_1 thus
%            Psi(y,lambda) = soft(y,lambda)
%            Important: if Psi is given, phi must also be given,
%                       so that the algorithm may also compute
%                       the objective function.
%
%  'StopCriterion' = type of stopping criterion to use
%                    0 = algorithm stops when the relative
%                        change in the number of non-zero
%                        components of the estimate falls
%                        below 'tolerance'
%                    1 = stop when the relative
%                       change in the objective function
%                       falls below 'tolerance'
%                    2 = stop when relative duality gap
%                       falls below 'tolerance'
%                    3 = stop when relative noise magnitude
%                       falls below 'tolerance'
%                    4 = stop when the objective function
%                        becomes equal or less than toleranceA.
%                    Default = 3
%
%  'Tolerance' = stopping threshold; Default = 0.01
%
%  'Debias'     = debiasing option: 1 = yes, 0 = no.
%                 Default = 0.
%
%  'ToleranceD' = stopping threshold for the debiasing phase:
%                 Default = 0.0001.
%                 If no debiasing takes place, this parameter,
%                 if present, is ignored.
%
%  'MaxiterA' = maximum number of iterations allowed in the
%               main phase of the algorithm.
%               Default = 1000
%
%  'MiniterA' = minimum number of iterations performed in the
%               main phase of the algorithm.
%               Default = 5
%
%  'MaxiterD' = maximum number of iterations allowed in the
%               debising phase of the algorithm.
%               Default = 200
%
%  'MiniterD' = minimum number of iterations to perform in the
%               debiasing phase of the algorithm.
%               Default = 5
%
%  'Initialization' must be one of {0,1,2,array}
%               0 -> Initialization at zero.
%               1 -> Random initialization.
%               2 -> initialization with A'*y.
%           array -> initialization provided by the user.
%               Default = 0;
%
%  'BB_variant' specifies which variant of Barzila-Borwein to use, or not.
%               0 -> don't use a BB rule - instead pick the starting alpha
%               based on the successful value at the previous iteration
%               1 -> standard BB choice  s'r/s's
%               2 -> inverse BB variant r'r/r's
%               Default = 1
%
%  'BB_cycle' specifies the cycle length  - the number of iterations between
%             recalculation of alpha. Requires integer value at least
%             1. Relevant only if a **nonmonotone BB rule** is used
%             (BB_variant = 1 or 2 and Monotone=0).
%             Default = 1
%
%  'Monotone' =  enforce monotonic decrease in f, or not?
%               any nonzero -> enforce monotonicity (overrides 'Safeguard')
%               0 -> don't enforce monotonicity.
%               Default = 0;
%
%  'Safeguard' = enforce a "sufficient decrease" over the largest
%               objective value of the past M iterations.
%               any nonzero -> safeguard
%               0 -> don't safeguard
%               Default = 0.
%
%  'M'        = number of steps to look back in the safeguarding process.
%               Ignored if Safeguard=0 or if Monotone is nonzero.
%               (positive integer. Default = 5)
%
%  'sigma'    = sigma value used in Safeguarding test for sufficient
%               decrease. Ignored unless 'Safeguard' is nonzero. Must be
%               in (0,1). Drfault: .01.
%
%  'Eta'      = factor by which alpha is multiplied within an iteration,
%               until a decrease in the objective function is
%               obtained.
%               Default = 2;
%
%  'Alpha_factor' = factor by which to reduce the successful value of
%                alpha at iteration k, to give the first value of alpha
%                to be tried at iteration k+1.
%                If a Barzilai-Borwein rule is specified (BB_variant > 0),
%                this parameter is ignored.
%                Default = 0.8;
%
%  'Continuation' = Continuation or not (1 or 0)
%                   Specifies the choice for a continuation scheme,
%                   in which we start with a large value of lambda, and
%                   then decrease lambda until the desired value is
%                   reached. At each value, the solution obtained
%                   with the previous values is used as initialization.
%                   Default = 0
%
% 'ContinuationSteps' = Number of steps in the continuation procedure;
%                       ignored if 'Continuation' equals zero.
%                       If -1, an adaptive continuation procedure is used.
%                       Default = -1.
%
% 'FirstlambdaFactor'  = Initial lambda value, if using continuation, is
%                     obtained by multiplying the given lambda by
%                     this factor. This parameter is ignored if
%                     'Continuation' equals zero or
%                     'ContinuationSteps' equals -1.
%                     Default = 10.
%
%  'True_x' = if the true underlying x is passed in
%                this argument, MSE plots are generated.
%
%  'AlphaMin' = the alphamin parameter of the BB method.
%               Default = 1e-30;
%
%  'AlphaMax' = the alphamax parameter of the BB method.
%               Default = 1e30;
%
%  'Verbose'  = work silently (0) or verbosely (1)
%
% ===================================================
% ============ Outputs ==============================
%   x = solution of the main algorithm
%
%   x_debias = solution after the debiasing phase;
%                  if no debiasing phase took place, this
%                  variable is empty, x_debias = [].
%
%
%   debias_start = iteration number at which the debiasing
%                  phase started. If no debiasing took place,
%                  this variable is returned as zero.
%
% ========================================================


% Set the defaults for the optional parameters
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;

STOPPING_DEFAULT = STOPPING_SPARSE_SUPPORT;
stoppingCriterion = STOPPING_DEFAULT;

tolerance = 0.01;
tolD = 0.0001;
debias = 0;
maxiter = 1000;
maxiter_debias = 200;
miniter = 5;
miniter_debias = 0;
bbVariant = 1;
bbCycle = 1;
enforceMonotone = 0;
enforceSafeguard = 0;
M = 5;
sigma = .01;
alphamin = 1e-30;
alphamax = 1e30;
verbose = 0;
continuation = 1;
cont_steps = 2;
psi_ok = 0;
xG = [];
% amount by which to increase alpha after an unsuccessful step
eta = 2.0;
% amount by which to decrease alpha between iterations, if a
% Barzilai-Borwein rule is not used to make the initial guess at each
% iteration.
alphaFactor = 0.8;

% Set the defaults for outputs that may not be computed
debias_start = 0;
x_debias = [];
e_debias = [];
x=[];
e=[];
iter = 1;
lambda = 1e-3;
isNonnegative = false;

% Read the optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'ISNONNEGATIVE'
                isNonnegative = varargin{i+1};
            case 'LAMBDA'
                lambda = varargin{i+1};
            case 'STOPPINGCRITERION'
                stoppingCriterion = varargin{i+1};
            case 'GROUNDTRUTH'
                xG = varargin{i+1};
            case 'TOLERANCE'
                tolerance = varargin{i+1};
            case 'TOLERANCED'
                tolD = varargin{i+1};
            case 'DEBIAS'
                debias = varargin{i+1};
            case 'MAXITERATION'
                maxiter = varargin{i+1};
            case 'MAXITERD'
                maxiter_debias = varargin{i+1};
            case 'MINITERA'
                miniter = varargin{i+1};
            case 'MINITERD'
                miniter_debias = varargin{i+1};
            case 'INITIALIZATION'
                x = varargin{i+1};
                e = y-A*x;
            case 'BB_VARIANT'
                bbVariant = varargin{i+1};
            case 'BB_CYCLE'
                bbCycle = varargin{i+1};
            case 'MONOTONE'
                enforceMonotone = varargin{i+1};
            case 'SAFEGUARD'
                enforceSafeguard = varargin{i+1};
            case 'M'
                M = varargin{i+1};
            case 'SIGMA'
                sigma = varargin{i+1};
            case 'ETA'
                eta = varargin{i+1};
            case 'ALPHA_FACTOR'
                alphaFactor = varargin{i+1};
            case 'CONTINUATION'
                continuation = varargin{i+1};
            case 'CONTINUATIONSTEPS'
                cont_steps = varargin{i+1};
            case 'ALPHAMIN'
                alphamin = varargin{i+1};
            case 'ALPHAMAX'
                alphamax = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

if stoppingCriterion==STOPPING_GROUND_TRUTH && isempty(xG)
    error('The stopping criterion must provide the ground truth value of x.');
end

maxiter_debias = min(maxiter_debias, ceil(maxiter/20));
%%%%%%%%%%%%%%

% it makes no sense to ask for a nonmonotone variant of a non-BB method
if ~enforceMonotone && bbVariant==0
    error('non-monotone, non-BBmethod requested');
end

[K,N] = size(A);
AT = A';

% Initialization
if isempty(x)
    x = zeros(N,1);
    e = zeros(K,1);
end

% if lambda is large enough, in the case of phi = l1, thus psi = soft,
% the optimal solution is the zero vector
aux = [AT*y; y];
max_lambda = max(abs(aux(:)));
firstlambdaFactor = 0.8*max_lambda / lambda;
if (lambda >= max_lambda) && (psi_ok==0)
    x = zeros(N,1);
    e = zeros(K,1);
    return
end

% define the indicator vector or matrix of nonzeros in x
nz_x = (abs([x; e]) > 10*eps);

% store given lambda, because we're going to change it in the
% continuation procedure
final_lambda = lambda;
% if we choose to use adaptive continuation, need to reset lambda to realmax to
% make things work (don't ask...)
if cont_steps == -1
    lambda = realmax;
end

% store given stopping criterion and threshold, because we're going
% to change them in the continuation procedure
final_tol = tolerance;

% set continuation factors

if (continuation && (cont_steps > 1))
    % If lambda is scalar, first check top see if the first factor is
    % too large (i.e., large enough to make the first
    % solution all zeros). If so, make it a little smaller than that.
    if numel(lambda) == 1
        if firstlambdaFactor*lambda >= max_lambda
            firstlambdaFactor = 0.5*max_lambda / lambda;
            if verbose
                fprintf(1,'\n setting parameter FirstlambdaFactor\n')
            end
        end
        cont_factors = 10.^[log10(firstlambdaFactor):...
            log10(1/firstlambdaFactor)/(cont_steps-1):0];
    end
else
    if ( ~continuation )
        cont_factors = 1;
        cont_steps = 1;
    end
end

keep_continuation = 1;
cont_loop = 1;

% loop for continuation
while keep_continuation
    
    % initialize the count of steps since last update of alpha
    % (for use in cyclic BB)
    iterThisCycle = 0;
    
    % Compute the initial residual and gradient
    resid =  A*x + e - y;
    gradqx = AT*resid; 
    gradqe = resid;
    
    if cont_steps == -1
        
        temp_lambda = max(final_lambda,0.2*max(abs([gradqx; gradqe])));
        
        if temp_lambda > lambda
            lambda = final_lambda;
        else
            lambda = temp_lambda;
        end
        
        if lambda == final_lambda
            currentStoppingCriterion = stoppingCriterion;
            tolerance = final_tol;
            keep_continuation = 0;
        else
            currentStoppingCriterion = STOPPING_OBJECTIVE_VALUE;
            tolerance = 1e-5;
        end
    else
        lambda = final_lambda * cont_factors(cont_loop);
        if cont_loop == cont_steps
            currentStoppingCriterion = stoppingCriterion;
            tolerance = final_tol;
            keep_continuation = 0;
        else
            currentStoppingCriterion = STOPPING_OBJECTIVE_VALUE;
            tolerance = 1e-5;
        end
    end
    
    if verbose
        fprintf('\n Regularization parameter lambda = %10.6e\n',lambda)
    end
    
    % compute and store initial value of the objective function
    % for this lambda
    alpha = 1.0;
    f = 0.5*(resid(:)'*resid(:)) + lambda * (norm(x,1) + norm(e,1));
    if enforceSafeguard
        f_lastM = f;
    end
    
    % initialization of alpha
    % alpha = 1/max(max(abs(du(:))),max(abs(dv(:))));
    % or just do a dumb initialization
    %alphas(iter) = alpha;
    
    % control variable for the outer loop and iteration counter
    keep_going = 1;
    
    while keep_going
        
        % compute gradient
        gradqx = AT*resid;
        gradqe = resid;
        
        % save current values
        prev_x = x;
        prev_e = e;
        prev_f = f;
        prev_resid = resid;
        
        % computation of step
        
        cont_inner = 1;
        while cont_inner
            if isNonnegative
                x = soft(max(prev_x - gradqx*(1/alpha),0), lambda/alpha);
                e = soft(max(prev_e - gradqe*(1/alpha),0), lambda/alpha);
            else
                x = soft(prev_x - gradqx*(1/alpha), lambda/alpha);
                e = soft(prev_e - gradqe*(1/alpha), lambda/alpha);
            end
            dx = x - prev_x;
            de = e - prev_e;
            Adx = A*dx;
            resid = prev_resid + Adx + de;
            f = 0.5*(resid(:)'*resid(:)) + lambda * norm(x,1);
            if enforceMonotone
                f_threshold = prev_f;
            elseif enforceSafeguard
                f_threshold = max(f_lastM) - 0.5*sigma*alpha*(dx'*dx + de'*de);
            else
                f_threshold = inf;
            end
            % f_threshold
            
            if f <= f_threshold
                cont_inner=0;
            else
                % not good enough, increase alpha and try again
                alpha = eta*alpha;
                if verbose
                    fprintf(1,' f=%10.6e, increasing alpha to %6.2e\n', f, alpha);
                end
            end
        end   % of while cont_inner
        
        if enforceSafeguard
            if length(f_lastM)<M+1
                f_lastM = [f_lastM f];
            else
                f_lastM = [f_lastM(2:M+1) f];
            end
        end
        
        % print stuff
        if verbose
            fprintf(1,'t=%4d, obj=%10.6e, alpha=%e  ', iter, f, alpha );
        end
        
        if bbVariant==1
            % standard BB choice of initial alpha for next step
            if iterThisCycle==0 || enforceMonotone==1
                dGd = Adx'*Adx + de'*de;
                alpha = min(alphamax,max(alphamin,dGd/(realmin+dx'*dx+de'*de)));
            end
        elseif bbVariant==2
            % alternative BB choice of initial alpha for next step
            if iterThisCycle==0 || enforceMonotone==1
                dGd = Adx'*Adx + de'*de;
                ATAdx= [AT*(Adx + de); Adx + de];
                dGGd = ATAdx'*ATAdx;
                alpha = min(alphamax,max(alphamin,dGGd/(realmin+dGd)));
            end
        else
            % reduce current alpha to get initial alpha for next step
            alpha = alpha * alphaFactor;
        end
        
        % update iteration counts, store results and times
        iter=iter+1;
        iterThisCycle=mod(iterThisCycle+1,bbCycle);
        % alphas(iter) = alpha;
        
        % compute stopping criteria and test for termination
        switch currentStoppingCriterion
            case STOPPING_GROUND_TRUTH
                keep_going = norm(xG-x)>tolerance;
            case STOPPING_SPARSE_SUPPORT
                % compute the stopping criterion based on the change
                % of the number of non-zero components of the estimate
                nz_x_prev = nz_x;
                nz_x = (abs([x; e])>10*eps);
                num_nz_x = sum(nz_x(:));
                num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
                if num_nz_x >= 1
                    criterionActiveSet = num_changes_active / num_nz_x;
                    keep_going = (criterionActiveSet > tolerance);
                end
            case STOPPING_OBJECTIVE_VALUE
                % compute the stopping criterion based on the relative
                % variation of the objective function.
                criterionObjective = abs(f-prev_f)/(prev_f);
                keep_going =  (criterionObjective > tolerance);
            case STOPPING_DUALITY_GAP
                % compute the "duality" stopping criterion - actually based on the
                % iterate PRIOR to the step just taken. Make it relative to the primal
                % function value.
                scaleFactor = norm([gradqx; gradqe],inf);
                w = lambda*resid / scaleFactor;
                criterionDuality = 0.5* (resid'*resid) + ...
                    lambda * norm(x,1) + 0.5*w(:)'*w(:) + y(:)'*w(:);
                % criterionDuality = criterionDuality / prev_f;
                keep_going = (criterionDuality > tolerance);
                if verbose
                    fprintf(1,'Duality = %e (target = %e)\n',...
                        criterionDuality , tolA)
                end
            case STOPPING_SUBGRADIENT
                error('Subgradient is not a valid stopping criterion for SpaRSA.');
            otherwise,
                error('Undefined stopping criterion.');
        end % end of the stopping criteria switch
        
        % overrule the stopping decision to ensure we take between miniter and
        % maxiter iterations
        if iter<=miniter
            % take no fewer than miniter...
            keep_going = 1;
        elseif iter > maxiter
            % and no more than maxiter iterations
            keep_going = 0;
        end
        
    end % end of the main loop of the GPBB algorithm (while keep_going)
    
    cont_loop = cont_loop + 1;
    
end % end of the continuation loop (while keep_continuation)

% If the 'Debias' option is set to 1, we try to
% remove the bias from the l1 penalty, by applying CG to the
% least-squares problem obtained by omitting the l1 term
% and fixing the zero coefficients at zero.

if (debias && (sum([x;e]~=0)~=0))
    if verbose
        fprintf(1,'\nStarting the debiasing phase...\n\n')
    end
    
    x_debias = x;
    e_debias = e;
    zeroindx = (x_debias~=0);
    zeroinde = (e_debias~=0);
    cont_debias_cg = 1;
    debias_start = iter;
    
    % calculate initial residual
    resid = A*x_debias + e_debias;
    resid = resid-y;
    prev_resid = eps*ones(size(resid));
    
    rvecx = (AT*resid).*zeroindx;
    rvece = resid.*zeroinde;
    
    % mask out the zeros
    rTr_cg = rvecx'*rvecx + rvece'*rvece;
    
    % set convergence threshold for the residual || RW x_debias - y ||_2
    tol_debias = tolD * rTr_cg;
    
    % initialize pvec
    pvecx = - rvecx;
    pvece = - rvece;
    
    % main loop
    while cont_debias_cg
        
        % calculate A*p = Wt * Rt * R * W * pvec
        RWpvec = A*pvecx + pvece;
        Apvecx = (AT*RWpvec).*zeroindx;
        Apvece = RWpvec.*zeroinde;
        
        % calculate alpha for CG
        alpha_cg = rTr_cg / (pvecx'* Apvecx+pvece'*Apvece);
        
        % take the step
        x_debias = x_debias + alpha_cg * pvecx;
        e_debias = e_debias + alpha_cg * pvece;
        resid = resid + alpha_cg * RWpvec;
        rvecx  = rvecx  + alpha_cg * Apvecx;
        rvece  = rvece  + alpha_cg * Apvece;
        
        rTr_cg_plus = rvecx'*rvecx + rvece'*rvece;
        beta_cg = rTr_cg_plus / rTr_cg;
        pvecx = -rvecx + beta_cg * pvecx;
        pvece = -rvece + beta_cg * pvece;
        
        rTr_cg = rTr_cg_plus;
        
        iter = iter+1;
        
        
        % in the debiasing CG phase, always use convergence criterion
        % based on the residual (this is standard for CG)
        if verbose
            fprintf(1,'t = %5d, debias resid = %13.8e, convergence = %8.3e\n', ...
                iter, resid(:)'*resid(:), rTr_cg / tol_debias);
        end
        cont_debias_cg = ...
            (iter-debias_start <= miniter_debias )| ...
            ((rTr_cg > tol_debias) & ...
            (iter-debias_start <= maxiter_debias));
        
    end
end

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = y./(y+T) .* x;
end