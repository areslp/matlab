function [x_t, nIter, timeSteps, errorSteps, nShrinkage] = SolveAMP(A, b, varargin)

% Solve
% min_x ||x||_1  s.t.  Ax = b
t0 = tic;

DEBUG = 0 ;
DISPLAY = 0 ;

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_INCREMENTS = 5 ;
STOPPING_FEASIBILITY = 6 ;
STOPPING_DEFAULT = STOPPING_FEASIBILITY;

stoppingCriterion = STOPPING_DEFAULT;

% Initialize parameters
[m,n] = size(A) ;
if m >= n
    error('The observation matrix must have more columns than rows.') ;
end

tol = 1e-6 ;

eps = 1e-6 ;

bNorm = norm(b) ;

maxIter = 10000 ;

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'initialization'
            disp(['Ignoring ' parameterName ' for this method.']) ;
            %             x0 = parameterValue;
            %             if ~all(size(x0)==[n,1])
            %                 error('The dimension of the initial x0 does not match.');
            %             end
        case 'groundtruth'
            xG = parameterValue;
        case 'mu'
            disp(['Ignoring ' parameterName ' for this method.']) ;
        case 'gamma'
            disp(['Ignoring ' parameterName ' for this method.']) ;
        case 'maxiteration'
            maxIter = parameterValue;
        case 'isnonnegative'
            disp(['Ignoring ' parameterName ' for this method.']) ;
        case 'tolerance'
            tol = parameterValue;
        case 'verbose'
            verbose = parameterValue;
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin
timeSteps = nan(1,maxIter) ;
errorSteps = nan(1,maxIter) ;

converged = 0 ;

nIter = 0 ;
nShrinkage = 0 ;

delta = m/n ;

tau_t = 0.1*norm(A'*b,inf) ;

deltaInv = n/m ;

x_t = zeros(n,1) ; z_t = b ;

f = norm(x_t,1) ;

nz_x = (abs(x_t)>eps*10);

while ~converged
    
    temp1 = A'*z_t ;
    
    if nIter == 0
        
        tau_t = 0.1*norm(temp1,inf) ;
        
    end
    
    temp2 = temp1 + x_t ;    
    
    x_tp1 = shrink(temp2, tau_t) ;
    
    residual = b - A*x_tp1 ;
    
    z_tp1 = residual + deltaInv*mean(shrink_der(temp2, tau_t))*z_t ;
    
    tau_tp1 = deltaInv*tau_t*mean(shrink_der(temp1 + x_tp1, tau_t)) ;
    
    nShrinkage = nShrinkage + 3 ;
    
    nIter = nIter + 1 ;
    
    timeSteps(nIter) = toc(t0) ;
    errorSteps(nIter) = norm(x_tp1-xG) ;
    
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            if norm(xG-x_tp1) < tol
                converged = 1 ;
            end
        case STOPPING_SUBGRADIENT
            error('Vanishing subgradient is not a valid stopping criterion for AMP.');
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            nz_x_prev = nz_x;
            nz_x = (abs(x_tp1)>eps*10);
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                converged = ~(criterionActiveSet > tol);
            end
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            prev_f = f;
            f = norm(x_tp1,1) ;
            criterionObjective = abs(f-prev_f)/(prev_f);
            converged =  ~(criterionObjective > tol);
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for AMP.');
        case STOPPING_INCREMENTS
            if norm(x_tp1 - x_t) < tol*norm(x_t)
                converged = 1 ;
            end            
        case STOPPING_FEASIBILITY
            if norm(residual) < tol*bNorm
                converged = 1 ;
            end
        case STOPPING_TIME
            converged = timeSteps(nIter) >= maxTime ;
        otherwise
            error('Undefined stopping criterion.');
    end
    
    if nIter >= maxIter && ~converged
        
        disp('Max. iterations reached.') ;
        converged = 1 ;
        
    end
    
    % if ~mod(nIter,2) && DISPLAY
    if DISPLAY > 1
        figure(100) ; clf ;
        stem(xtp1) ;
        title(['Iteration ' num2str(nIter)]) ;
        pause(0.1) ;
        
    end
    
    if DISPLAY > 0
        
        disp(['Iteration ' num2str(nIter) ' ||x||_0 ' num2str(sum(double(abs(x_tp1)>0)))]) ;
        
    end
    
    x_t = x_tp1 ;
    z_t = z_tp1 ;
    tau_t = tau_tp1 ;
    
end
timeSteps = timeSteps(1:nIter) ;
errorSteps = errorSteps(1:nIter) ;

function Y = shrink(X, alpha)
    
Y = sign(X).*max(abs(X)-alpha,0) ;

function Y = shrink_der(X, alpha)

Y = double(abs(X) > alpha) ;