% The following primal-dual interior-point algorithm is modified from l1eq_pd.m
%
% Solve
% min_x ||x||_1  s.t.  Ax = b
%
% Recast as linear program
% min_{x,u} sum(u)  s.t.  -u <= x <= u,  Ax=b
% and use primal-dual interior point method
%
% Usage: xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - A KxN matrix.
%
%
% b - Kx1 vector of observations.
%
% pdtol - Tolerance for primal-dual algorithm (algorithm terminates if
%     the duality gap is less than pdtol).
%     Default = 1e-3.
%
% pdmaxiter - Maximum number of primal-dual iterations.
%     Default = 50.
%
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [x_out, e_out, pditer] = SolvePDIPA_CBM(A, b, varargin)

DEBUG = 1;

STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_DUALITY_GAP;

stoppingCriterion = STOPPING_DEFAULT;

tolerance = 1e-3;
pdmaxiter = 50;
[K, N] = size(A);
n = K+N;
x0 = zeros(n,1);
At = A';

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
            x0 = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'lambda'
            lambda = parameterValue;
        case 'maxiteration'
            pdmaxiter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'tolerance'
            tolerance = parameterValue;
        case 'verbose'
            verbose = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

alpha = 0.01;
beta = 0.5;
mu = 10;

gradf0 = [zeros(n,1); ones(n,1)];

x = x0;
u = 1.01*max(abs(x))*ones(n,1) + 1e-2;

fu1 = x - u;
fu2 = -x - u;

lamu1 = -1./fu1;
lamu2 = -1./fu2;
v = -A*(lamu1(1:N)-lamu2(1:N)) - (lamu1(N+1:end)-lamu2(N+1:end));
Atv = [At*v; v];
rpri = A*x(1:N) + x(N+1:end) - b;

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*N/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(n,1)];
resnorm = norm([rdual; rcent; rpri]);

pditer = 0;
while (pditer < pdmaxiter)
    
    pditer = pditer + 1;
    
    w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;
    w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
    w3 = -rpri;
    
    sig1 = -lamu1./fu1 - lamu2./fu2;
    sig2 = lamu1./fu1 - lamu2./fu2;
    sigx = sig1 - sig2.^2./sig1;
    
    if any(sigx==0)
        sigx = sigx + 100*eps;
    end
    
    H11p = -A*diag(1./sigx(1:N))*At-diag(1./sigx(N+1:end));
    w1p = (w1./sigx - w2.*sig2./(sigx.*sig1));
    w1p = w3 - A*w1p(1:N) - w1p(N+1:end);
    [dv,hcond] = linsolve(H11p,w1p);
    if (hcond < 1e-14)
        if DEBUG>0
            disp('Primal-dual: Matrix ill-conditioned.  Returning previous iterate.');
        end
        x_out = x(1:N);
        e_out = x(N+1:end);
        return
    end
    dx = (w1 - w2.*sig2./sig1 - [At*dv; dv])./sigx;
    Adx = A*dx(1:N)+dx(N+1:end);
    Atdv = [At*dv; dv];
    
    du = (w2 - sig2.*dx)./sig1;
    
    dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;
    dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2;
    
    % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
    indp = find(dlamu1 < 0);  indn = find(dlamu2 < 0);
    s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);
    indp = find((dx-du) > 0);  indn = find((-dx-du) > 0);
    s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);
    
    % backtracking line search
    backiter = 0;
    xp = x + s*dx;  up = u + s*du;
    vp = v + s*dv;  Atvp = Atv + s*Atdv;
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
    fu1p = xp - up;  fu2p = -xp - up;
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(n,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    while(norm([rdp; rcp; rpp]) > (1-alpha*s)*resnorm)
        s = beta*s;
        xp = x + s*dx;  up = u + s*du;
        vp = v + s*dv;  Atvp = Atv + s*Atdv;
        lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
        fu1p = xp - up;  fu2p = -xp - up;
        rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(n,1)];
        rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
        rpp = rpri + s*Adx;
        backiter = backiter+1;
        if (backiter > 32)
            if DEBUG>0
                disp('Stuck backtracking, returning last iterate.')
            end
            x_out = x(1:N);
            e_out = x(N+1:end);
            return
        end
    end
    
    v = vp;  Atv = Atvp;
    lamu1 = lamu1p;  lamu2 = lamu2p;
    fu1 = fu1p;  fu2 = fu2p;
    
    % surrogate duality gap
    sdg = -(fu1'*lamu1 + fu2'*lamu2);
    tau = mu*2*n/sdg;
    rpri = rpp;
    rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
    rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(n,1)];
    resnorm = norm([rdual; rcent; rpri]);
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            done = norm(xp(1:N)-xG)<tolerance;
        case STOPPING_SPARSE_SUPPORT
            error('Sparse support is not a valid stopping criterion for BP.');
        case STOPPING_DUALITY_GAP
            done = (sdg < tolerance);
        case STOPPING_OBJECTIVE_VALUE
            % continue if not yeat reached target value tolA
            error('Objective value is not a valid stopping criterion for BP.');
        case STOPPING_SUBGRADIENT
            error('Subgradient is not a valid stopping criterion for BP.');
        otherwise,
            error('Undefined stopping criterion');
    end % end of the stopping criteria switch
    
    if done || norm(x-xp)<100*eps
        break;
    end
    
    % next iteration
    x = xp;  u = up;
end

x_out = xp(1:N);
e_out = xp(N+1:end);