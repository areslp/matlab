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

%% This function is modified from Matlab code proximal_gradient_bp

function [x_hat,nIter, timeSteps, errorSteps] = SolveFISTA(A,b, varargin)

% b - m x 1 vector of observations/data (required input)
% A - m x n measurement matrix (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maxilambdam number of iterations
%         - DEFAULT 10000, if omitted or -1.
% lineSearchFlag - 1 if line search is to be done every iteration
%                - DEFAULT 0, if omitted or -1.
% continuationFlag - 1 if a continuation is to be done on the parameter lambda
%                  - DEFAULT 1, if omitted or -1.
% eta - line search parameter, should be in (0,1)
%     - ignored if lineSearchFlag is 0.
%     - DEFAULT 0.9, if omitted or -1.
% lambda - relaxation parameter
%    - ignored if continuationFlag is 1.
%    - DEFAULT 1e-3, if omitted or -1.
% outputFileName - Details of each iteration are dumped here, if provided.
%
% x_hat - estimate of coeeficient vector
% numIter - number of iterations until convergence
%
%
% References
% "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
% "An Accelerated Proximal Gradient Algorithm for Nuclear Norm Regularized Least Squares problems", K.-C. Toh and S. Yun, preprint 2009.
%
% Arvind Ganesh, Summer 2009. Questions? abalasu2@illinois.edu

t0 = tic ;

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_SUBGRADIENT;

stoppingCriterion = STOPPING_DEFAULT;
maxIter = 10000 ;
tolerance = 1e-3;
[m,n] = size(A) ;
x0 = zeros(n,1) ;
xG = [];

%% Initializing optimization variables
t_k = 1 ; 
t_km1 = 1 ;
L0 = 1 ;
G = A'*A ;
nIter = 0 ;
c = A'*b ;
lambda0 = 0.5*L0*norm(c,inf) ;
eta = 0.95 ;
% lambda_bar = 1e-10*lambda0 ;
lambda_bar = 1e-6;
% xk = zeros(n,1) ;
xk = A\b ;
lambda = lambda0 ;
L = L0 ;
beta = 1.5 ;

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
            tolerance = parameterValue;
        case 'linesearchflag'
            lineSearchFlag = parameterValue;
        case 'lambda'
            lambda_bar = parameterValue;
        case 'maxiteration'
            maxIter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'continuationflag'
            continuationFlag = parameterValue;
        case 'initialization'
            xk = parameterValue;
            if ~all(size(xk)==[n,1])
                error('The dimension of the initial xk does not match.');
            end
        case 'eta'
            eta = parameterValue;
            if ( eta <= 0 || eta >= 1 )
                disp('Line search parameter out of bounds, switching to default 0.9') ;
                eta = 0.9 ;
            end
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = nan(1,maxIter) ;
errorSteps = nan(1,maxIter) ;

if stoppingCriterion==STOPPING_GROUND_TRUTH && isempty(xG)
    error('The stopping criterion must provide the ground truth value of x.');
end

keep_going = 1 ;
nz_x = (abs(xk)> eps*10);
f = 0.5*norm(b-A*xk)^2 + lambda_bar * norm(xk,1);
xkm1 = xk;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;
    
    yk = xk + ((t_km1-1)/t_k)*(xk-xkm1) ;
    
    stop_backtrack = 0 ;
    
    temp = G*yk - c ; % gradient of f at yk
    
    while ~stop_backtrack
        
        gk = yk - (1/L)*temp ;
        
        xkp1 = soft(gk,lambda/L) ;
        
        temp1 = 0.5*norm(b-A*xkp1)^2 ;
        temp2 = 0.5*norm(b-A*yk)^2 + (xkp1-yk)'*temp + (L/2)*norm(xkp1-yk)^2 ;
        
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*beta ;
        end
        
    end
    
    % disp(['Iter ' num2str(nIter) ' ||x||_0 ' num2str(sum(abs(xkp1) > 0)) ' err ' num2str(norm(xG-xkp1))]) ;
    timeSteps(nIter) = toc(t0) ;
    errorSteps(nIter) = norm(xkp1-xG) ;
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            keep_going = norm(xG-xkp1)>tolerance;
        case STOPPING_SUBGRADIENT
            sk = L*(yk-xkp1) + G*(xkp1-yk) ;
            keep_going = norm(sk) > tolerance*L*max(1,norm(xkp1));
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            nz_x_prev = nz_x;
            nz_x = (abs(xkp1)>eps*10);
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                keep_going = (criterionActiveSet > tolerance);
            end
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            prev_f = f;
            f = 0.5*norm(b-A*xkp1)^2 + lambda_bar * norm(xk,1);
            criterionObjective = abs(f-prev_f)/(prev_f);
            keep_going =  (criterionObjective > tolerance);
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for PGBP.');
        case STOPPING_TIME
            keep_going = timeSteps(nIter) < maxTime ;
        otherwise
            error('Undefined stopping criterion.');
    end
    
    lambda = max(eta*lambda,lambda_bar) ;
    
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    t_km1 = t_k ;
    t_k = t_kp1 ;
    xkm1 = xk ;
    xk = xkp1 ;
end

x_hat = xk ;
timeSteps = timeSteps(1:nIter) ;
errorSteps = errorSteps(1:nIter) ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end