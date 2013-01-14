function [x, e, nIter, timeSteps, errorSteps, idSteps] = SolveDALM_CBM(A, b, varargin)

t0 = tic ;

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_INCREMENTS = 5 ;
STOPPING_DEFAULT = STOPPING_DUALITY_GAP;

stoppingCriterion = STOPPING_DEFAULT;

maxIter = 500;
nu = 1 ;
VERBOSE = 0 ;
tol = 1e-5 ;

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
        case 'nu'
            nu = parameterValue ;
        case 'maxtime'
            maxTime = parameterValue;
        case 'recdata'
            recData = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

% Initialize parameters
[m,n] = size(A) ;

Av = [A nu * eye(m)] / sqrt(1 + nu^2);
bv = 1 / sqrt(1+nu^2) * b;

G = (A * A' + nu^2 * eye(m))/ (1+nu^2);
invG_Av = G \ Av;
invG_bv = G \ bv;

beta = norm(bv,1)/m;
betaInv = 1/beta ;

nIter = 0 ;

if VERBOSE
    disp(['beta is: ' num2str(beta)]);
end

y = zeros(m,1);
x = zeros(n,1);     e = zeros(m,1);
xv = [x; e];
z = zeros(m+n,1);

converged_main = 0 ;

timeSteps = nan(1,maxIter) ;
errorSteps = nan(1,maxIter) ;
idSteps = nan(1,maxIter);

temp = Av' * y;
while ~converged_main
    
    nIter = nIter + 1 ;
    
    x_old = x;
    e_old = e;
    
    %update z
    temp1 = temp + xv * betaInv;
    z = sign(temp1) .* min(1,abs(temp1));
    
    %compute Av' * y
    y = invG_Av * (z - xv * betaInv) + invG_bv * betaInv;
    %temp = Av_invG_Av * (z - xv * betaInv) + Av_invG_bv * betaInv;
    %temp = Av' * (invG * (Av * (z - xv * betaInv))) + Av_invG_bv *
    %betaInv;
    temp = Av' * y;
    
    %update x_v
    xv = xv - beta * (z - temp);
    
    x = xv(1:n);
    e = xv(n+1:end);
    
    if VERBOSE && mod(nIter, 50) == 0
        
        disp(['Iteration ' num2str(nIter)]) ;
        disp([norm(x-x_old)/norm(x_old) norm(e-e_old)/norm(e_old)]);
        
        figure(1);
        subplot(3,1,1);
        plot(x);
        title('x');
        subplot(3,1,2);
        plot(e);
        title('e');
        subplot(3,1,3);
        plot(z);
        title('z');
        pause;
        
    end
    
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
            if abs(norm(xv,1)- y.'*b)<tol
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
    
    if ~converged_main && nIter >= maxIter
        disp('Maximum Iterations Reached') ;
        converged_main = 1 ;
    end
    
end

x = xv(1:n);
e = xv(n+1:end);

timeSteps = timeSteps(1:nIter) ;
errorSteps = errorSteps(1:nIter) ;
idSteps = idSteps(1:nIter);