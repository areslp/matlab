% Copyright ©2010. The Regents of the University of California (Regents). 
% All Rights Reserved. Contact The Office of Technology Licensing, 
% UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, 
% (510) 643-7201, for commercial licensing opportunities.

% Authors: Arvind Ganesh, Allen Y. Yang, Zihan Zhou.
% Contact: Allen Y. Yang, Department of EECS, University of California,
% Berkeley. <yang@eecs.berkeley.edu>

% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, 
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, 
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF 
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, 
% PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO 
% PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

%% The following primal-dual interior-point algorithm is modified from l1eq_pd.m
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
% A - Either a handle to a function that takes a N vector and returns a K
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
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
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [xp, pditer, timeSteps, errorSteps] = SolvePDIPA(A, b, varargin)

t0 = tic ;
DEBUG = 1;

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_DUALITY_GAP;

stoppingCriterion = STOPPING_DEFAULT;

tolerance = 1e-3;
pdmaxiter = 50;
cgtol = 1e-8;
cgmaxiter = 200;
N = size(A,2);
x0 = zeros(N,1);
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
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = nan(1,pdmaxiter) ;
errorSteps = nan(1,pdmaxiter) ;

alpha = 0.01;
beta = 0.5;
mu = 10;

gradf0 = [zeros(N,1); ones(N,1)];

x = x0;
u = 1.01*max(abs(x))*ones(N,1) + 1e-2;

fu1 = x - u;
fu2 = -x - u;

lamu1 = -1./fu1;
lamu2 = -1./fu2;
v = -A*(lamu1-lamu2);
Atv = At*v;
rpri = A*x - b;

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*N/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
resnorm = norm([rdual; rcent; rpri]);

pditer = 0;
f = norm(x,1);
while (pditer < pdmaxiter)
    
    pditer = pditer + 1;
    
    w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;
    w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
    w3 = -rpri;
    
    sig1 = -lamu1./fu1 - lamu2./fu2;
    sig2 = lamu1./fu1 - lamu2./fu2;
    sigx = sig1 - sig2.^2./sig1;
    
    if any(sigx==0)
        sigx = sigx + 1e-10;
    end
    
    H11p = -A*diag(1./sigx)*At;
    w1p = w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1));
    [dv,hcond] = linsolve(H11p,w1p);
    if (hcond < 1e-14)
        if DEBUG>0
            disp('Primal-dual: Matrix ill-conditioned.  Returning previous iterate.');
        end
        xp = x;
        timeSteps = timeSteps(1:pditer-1) ;
        errorSteps = errorSteps(1:pditer-1) ;
        return;
    end
    dx = (w1 - w2.*sig2./sig1 - At*dv)./sigx;
    Adx = A*dx;
    Atdv = At*dv;
    
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
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    while(norm([rdp; rcp; rpp]) > (1-alpha*s)*resnorm)
        s = beta*s;
        xp = x + s*dx;  up = u + s*du;
        vp = v + s*dv;  Atvp = Atv + s*Atdv;
        lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
        fu1p = xp - up;  fu2p = -xp - up;
        rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
        rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
        rpp = rpri + s*Adx;
        backiter = backiter+1;
        if (backiter > 32)
            if DEBUG>0
                disp('Stuck backtracking, returning last iterate.')
            end
            xp = x;
            timeSteps = timeSteps(1:pditer-1) ;
            errorSteps = errorSteps(1:pditer-1) ;
            return;
        end
    end
    
    v = vp;  Atv = Atvp;
    lamu1 = lamu1p;  lamu2 = lamu2p;
    fu1 = fu1p;  fu2 = fu2p;
    
    % surrogate duality gap
    sdg = -(fu1'*lamu1 + fu2'*lamu2);
    tau = mu*2*N/sdg;
    rpri = rpp;
    rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
    rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
    resnorm = norm([rdual; rcent; rpri]);
    
    timeSteps(pditer) = toc(t0) ;
    errorSteps(pditer) = norm(xp-xG) ;
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            done = norm(xp-xG)<tolerance;
        case STOPPING_SPARSE_SUPPORT
            error('Sparse support is not a valid stopping criterion for BP.');
        case STOPPING_DUALITY_GAP
            done = (sdg < tolerance);
        case STOPPING_OBJECTIVE_VALUE
            prev_f = f;
            f = norm(x,1);
            criterionObjective = abs(f-prev_f)/(prev_f);
            done = (criterionObjective<tolerance);
        case STOPPING_SUBGRADIENT
            error('Subgradient is not a valid stopping criterion for BP.');
        case STOPPING_TIME
            done = timeSteps(pditer) >= maxTime ;
        otherwise,
            error('Undefined stopping criterion');
    end % end of the stopping criteria switch
    
    if done || norm(x-xp)<100*eps
        break;
    end
    
    % next iteration
    x = xp;  u = up;
end

timeSteps = timeSteps(1:pditer) ;
errorSteps = errorSteps(1:pditer) ;
