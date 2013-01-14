%% This function is modified from Matlab Package l1_ls

function [x_out, e_out, ntiter ,timeSteps, errorSteps, idSteps, status] = SolveL1LS_CBM(A,y,varargin)
%
% l1-Regularized Least Squares Problem Solver
%
%   l1_ls solves problems of the following form:
%
%       minimize ||A*x-y||^2 + lambda*sum|x_i|,
%
%   where A and y are problem data and x is variable (described below).
%
%
% INPUT
%   A       : mxn matrix; input data. columns correspond to features.
%
%   At      : nxm matrix; transpose of A.
%   m       : number of examples (rows) of A
%   n       : number of features (column)s of A
%
%   y       : m vector; outcome.
%   lambda  : positive scalar; regularization parameter
%
%   tar_gap : relative target duality gap (default: 1e-3)
%   quiet   : boolean; suppress printing message when true (default: false)
%
%   (advanced arguments)
%       eta     : scalar; parameter for PCG termination (default: 1e-3)
%       pcgmaxi : scalar; number of maximum PCG iterations (default: 5000)
%
% OUTPUT
%   x       : n vector; classifier
%   status  : string; 'Solved' or 'Failed'
%
%
% USAGE EXAMPLES
%   [x,status] = l1_ls(A,y,lambda);
%   [x,status] = l1_ls(A,At,m,n,y,lambda,0.001);
%

% AUTHOR    Kwangmoo Koh <deneb1@stanford.edu>
% UPDATE    Apr 8 2007
%
% COPYRIGHT 2008 Kwangmoo Koh, Seung-Jean Kim, and Stephen Boyd

%------------------------------------------------------------
%       INITIALIZE
%------------------------------------------------------------
t0 = tic ;

% IPM PARAMETERS
MU              = 2;        % updating parameter of t
MAX_NT_ITER     = 400;      % maximum IPM (Newton) iteration

% LINE SEARCH PARAMETERS
ALPHA           = 0.01;     % minimum fraction of decrease in the objective
BETA            = 0.5;      % stepsize decrease factor

At = A';
[m,N] = size(A);
n = m + N;
reltol = 1e-3;
quiet = true;
isNonnegative = false;

lambda = 1e-1;
lambda_bar = 1e-6;
rho = .1;

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_DUALITY_GAP;

stoppingCriterion = STOPPING_DEFAULT;

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'lambda'
            lambda = parameterValue;
        case 'maxiteration'
            MAX_NT_ITER = parameterValue;
        case 'tolerance'
            reltol = parameterValue;
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'maxtime'
            maxTime = parameterValue;
        case 'recdata'
            recData = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

if isNonnegative == true
    error('The current implementation does not support explicit nonnegativity contraint.');
end
    
MAX_LS_ITER     = ceil(MAX_NT_ITER/4);      % maximum backtracking line search iteration

% VARIABLE ARGUMENT HANDLING
t = min(max(1,1/lambda),2*n/1e-3);
eta = 1e-3;
pcgmaxi = 5000;
x = zeros(n,1);
u = ones(n,1);
objective = lambda*norm(x,1) + 1/2*norm(y-A*x(1:N)-x(N+1:end))^2;
f = [x-u;-x-u];

timeSteps = nan(1,MAX_NT_ITER+1) ;
errorSteps = nan(1,MAX_NT_ITER+1) ;
idSteps = nan(1,MAX_NT_ITER+1) ;

% RESULT VARIABLES
s   = Inf; pitr  = 0 ; dobj  =-Inf;

ntiter  = 0; lsiter  = 0; dxu =  zeros(2*n,1);

% diagxtx = diag(At*A);
diagxtx = 2*ones(n,1);

%------------------------------------------------------------
%               MAIN LOOP
%------------------------------------------------------------
previous_gap = 10;
oldX = x;
for ntiter = 0:MAX_NT_ITER
    
    z = A*x(1:N)+x(N+1:end)-y;
    
    %------------------------------------------------------------
    %       CALCULATE DUALITY GAP
    %------------------------------------------------------------
    
    nu = 2*z;
    
    maxAnu = norm([At*nu; nu],inf);
    if (maxAnu > lambda)
        nu = nu*lambda/maxAnu;
    end
    pobj  =  0.5*(z'*z)+lambda*norm(x,1);
    dobj  =  max(-0.5*(nu'*nu)-nu'*y,dobj);
    gap   =  pobj - dobj;
    
    %------------------------------------------------------------
    %       UPDATE t
    %------------------------------------------------------------
    if (s >= 0.5)
        t = max(min(2*n*MU/gap, MU*t), t);
    end
    
    %------------------------------------------------------------
    %       CALCULATE NEWTON STEP
    %------------------------------------------------------------
    
    q1 = 1./(u+x);          q2 = 1./(u-x);
    d1 = (q1.^2+q2.^2)/t;   d2 = (q1.^2-q2.^2)/t;
    
    
    % calculate gradient
    gradphi = [2*[At*z; z]-(q1-q2)/t; lambda*ones(n,1)-(q1+q2)/t];
    
    % calculate vectors to be used in the preconditioner
    prb     = diagxtx+d1;
    prs     = prb.*d1-(d2.^2);
    
    % set pcg tolerance (relative)
    normg   = norm(gradphi);
    pcgtol  = min(1e-1,eta*gap/min(1,normg));
    
    if (ntiter ~= 0 && pitr == 0) pcgtol = pcgtol*0.1; end
    
    [dxu,pflg,prelres,pitr,presvec] = ...
        pcg(@AXfunc_l1_ls,-gradphi,pcgtol,pcgmaxi,@Mfunc_l1_ls,...
        [],dxu,A,At,d1,d2,d1./prs,d2./prs,prb./prs);
    
    if (pflg == 1) pitr = pcgmaxi; end
    
    dx  = dxu(1:n);
    du  = dxu(n+1:end);
    
    %------------------------------------------------------------
    %   BACKTRACKING LINE SEARCH
    %------------------------------------------------------------
    phi = z'*z+lambda*sum(u)-sum(log(-f))/t;
    s = 1.0;
    gdx = gradphi'*dxu;
    for lsiter = 1:MAX_LS_ITER
        newx = x+s*dx; newu = u+s*du;
        newf = [newx-newu;-newx-newu];
        if (max(newf) < 0)
            newz   =  A*newx(1:N) + newx(N+1:end) - y;
            newphi =  newz'*newz+lambda*sum(newu)-sum(log(-newf))/t;
            if (newphi-phi <= ALPHA*s*gdx)
                break;
            end
        end
        s = BETA*s;
    end
    if (lsiter == MAX_LS_ITER) break; end % exit by BLS
    
    x = newx; u = newu; f = newf;
    
    timeSteps(ntiter+1) = toc(t0) ;    
    errorSteps(ntiter+1) = norm(x(1:N)-xG) ;
    sci = zeros(1,recData.nSubj);
    norm_xk = norm(x(1:N),1);
    for i = 1:recData.nSubj
        sci(i) = (recData.nSubj * norm(x(((i-1)*recData.nImg+1):i*recData.nImg), 1)/norm_xk - 1) ...
            / (recData.nSubj - 1);
    end
    [dontCare, curId] = max(sci);
    idSteps(ntiter+1) = (curId == recData.id);
    
    %------------------------------------------------------------
    %   STOPPING CRITERION
    %------------------------------------------------------------
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            keep_going = norm(xG-x(1:N))>reltol;
        case STOPPING_OBJECTIVE_VALUE
            prev_objective = objective;
            objective = lambda*norm(x,1) + 1/2*norm(y-A*x(1:N)-x(N+1:end))^2;
            keep_going = abs((objective-prev_objective)/prev_objective) > reltol;
        case STOPPING_DUALITY_GAP
            keep_going = (abs(gap - previous_gap)/previous_gap >= reltol);
            previous_gap = gap;
        case STOPPING_SPARSE_SUPPORT
            error('Sparse support is not a valid stopping criterion for L1LS.');
        case STOPPING_SUBGRADIENT
            error('Subgradient is not a valid stopping criterion for L1LS.');
        case STOPPING_TIME
            keep_going = timeSteps(ntiter+1) < maxTime ;
        otherwise
            error('Undefined stopping criterion');
    end
    
    if ~keep_going
        status  = 'Solved';
        break;
    end
    
    if norm(oldX-x)<1e-6 * norm(oldX)
        lambda = max(lambda * rho, lambda_bar);
        % RESULT VARIABLES
        dobj  =-Inf;
    end
    oldX = x;

end


%------------------------------------------------------------
%       ABNORMAL TERMINATION (FALL THROUGH)
%------------------------------------------------------------
if keep_going
    status = 'Failed';
    
    if (lsiter == MAX_LS_ITER)
        % failed in backtracking linesearch.
        if (~quiet) disp('MAX_LS_ITER exceeded in BLS'); end
    elseif (ntiter == MAX_NT_ITER)
        % fail to find the solution within MAX_NT_ITER
        if (~quiet) disp('MAX_NT_ITER exceeded.'); end
    elseif (~quiet) 
        disp('Search stuck at local minimum.');
    end
end

x_out = x(1:N);
e_out = x(N+1:end);

timeSteps = timeSteps(1:ntiter+1) ;
errorSteps = errorSteps(1:ntiter+1) ;
idSteps = idSteps(1:ntiter+1) ;

%------------------------------------------------------------
%       COMPUTE AX (PCG)
%------------------------------------------------------------
function [y] = AXfunc_l1_ls(x,A,At,d1,d2,p1,p2,p3)
%
% y = hessphi*[x1;x2],
%
% where hessphi = [A'*A*2+D1 , D2;
%                  D2        , D1];

n  = length(x)/2;
[m,N] = size(A);
x1 = x(1:n);
x2 = x(n+1:end);

y = (A*x1(1:N)+x1(N+1:end))*2;
y = [[At*y; y]+d1.*x1+d2.*x2; d2.*x1+d1.*x2];

%------------------------------------------------------------
%       COMPUTE P^{-1}X (PCG)
%------------------------------------------------------------
function [y] = Mfunc_l1_ls(x,A,At,d1,d2,p1,p2,p3)
%
% y = P^{-1}*x,
%

n  = length(x)/2;
x1 = x(1:n);
x2 = x(n+1:end);

y = [ p1.*x1-p2.*x2;...
    -p2.*x1+p3.*x2];

