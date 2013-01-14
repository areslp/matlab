function [x, e, nIter, timeSteps, errorSteps, idSteps] = SolvePALM_CBM(A, b, varargin)

t0 = tic ;
DEBUG = 0 ;

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_INCREMENTS = 5 ;
STOPPING_DEFAULT = STOPPING_DUALITY_GAP;

stoppingCriterion = STOPPING_DEFAULT;

% Initialize parameters
[m,n] = size(A) ;

tol = 5e-2 ;
tol_int = 1e-6 ;

G = A'*A ;
opts.disp = 0;
tau = eigs(G,1,'lm',opts);
tauInv = 1/tau ;

nIter = 0 ;

mu = 2 *m / norm(b,1);

lambda = zeros(m,1);
x = zeros(n,1) ;
e = b ;

converged_main = 0 ;

maxIter = 200 ;
maxIter_apg = 400;

nz_x = (abs([x; e])>eps*10);
f = norm([x;e],1) ;

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' lambdast be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'tolerance'
            tol = parameterValue;
        case 'maxiteration'
            maxIter = parameterValue;
        case 'maxtime'
            maxTime = parameterValue;
        case 'recdata'
            recData = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = nan(1,maxIter) ;
errorSteps = nan(1,maxIter) ;
idSteps = nan(1,maxIter) ;

while ~converged_main
    
    muInv = 1/mu ;
    lambdaScaled = muInv*lambda ;
    
    nIter = nIter + 1 ;
    
    e_old_main = e ;
    x_old_main = x ;
    
    temp2 = b + lambdaScaled ;
    temp = temp2 - A*x ;
    
    e = shrink(temp,muInv) ;
    
    converged_apg = 0 ;
    
    temp1 = A'*(e - temp2) ;
    
    nIter_apg = 0 ;
    
    t1 = 1 ; z = x ;
    
    muTauInv = muInv*tauInv ;
    
    Gx = G * x;
    Gz = Gx;
    
    while ~converged_apg
        
        nIter_apg = nIter_apg + 1 ;
        
        x_old_apg = x ;
        Gx_old = Gx;
        
        temp = z - tauInv*(temp1 + Gz) ;
        
        %x = shrink(temp, muTauInv) ;
        x = sign(temp) .* max(abs(temp)-muTauInv, 0);
        
        Gx = G * x;
        
        s = tau * (z - x) + Gx - Gz;
        if norm(s) < tol_int * tau * max(1,norm(x))
            converged_apg = 1;
        end
        
        if nIter_apg >= maxIter_apg
            converged_apg = 1 ;
        end
        
        t2 = (1+sqrt(1+4*t1*t1))/2 ;
        z = x + ((t1-1)/t2)*(x-x_old_apg) ;
        Gz = Gx + ((t1-1)/t2) * (Gx - Gx_old);
        t1 = t2 ;
        
    end
    
    lambda = lambda + mu*(b - A*x - e) ;
    
    timeSteps(nIter) = toc(t0) ;   
    errorSteps(nIter) = norm(x-xG) ;
    sci = zeros(1,recData.nSubj);
    norm_xk = norm(x,1);
    for i = 1:recData.nSubj
        sci(i) = (recData.nSubj * norm(x(((i-1)*recData.nImg+1):i*recData.nImg), 1)/norm_xk - 1) ...
            / (recData.nSubj - 1);
    end
    [dontCare, curId] = max(sci);
    idSteps(nIter) = (curId == recData.id);
        
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            if norm(xG-x) < tol
                converged_main = 1 ;
            end
        case STOPPING_SUBGRADIENT
            error('Duality gap is not a valid stopping criterion for ALM.');
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            nz_x_prev = nz_x;
            nz_x = (abs([x; e])>eps*10);
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                converged_main = ~(criterionActiveSet > tol);
            end
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            prev_f = f;
            f = norm([x ; e],1);
            criterionObjective = abs(f-prev_f)/(prev_f);
            converged_main =  ~(criterionObjective > tol);
        case STOPPING_DUALITY_GAP
            if abs(norm(x,1)- lambda.'*b)<tol
                converged_main = 1;
            end
        case STOPPING_INCREMENTS
            if norm([x_old_main ; e_old_main] - [x ; e]) < tol*norm([x_old_main ; e_old_main])
                converged_main = 1 ;
            end
        case STOPPING_TIME
            converged_main = timeSteps(nIter) >= maxTime ;
        otherwise
            error('Undefined stopping criterion.');
    end    
    
%     if ~converged_main && norm(x_old_main-x) < 100*eps
%         if DEBUG
%             disp('The iteration is stuck.') ;
%         end
%         converged_main = 1 ;
%         
%     end

    
    if ~converged_main && nIter >= maxIter
        if DEBUG
            disp('Maximum Iterations Reached') ;
        end
        converged_main = 1 ;
        
    end
        
end
timeSteps = timeSteps(1:nIter) ;
errorSteps = errorSteps(1:nIter) ;

function Y = shrink(X, alpha)
    
Y = sign(X).*max(abs(X)-alpha,0);