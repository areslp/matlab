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

function [x, nIter, timeSteps, errorSteps] = SolvePALM(A, b, varargin)

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
STOPPING_DEFAULT = STOPPING_INCREMENTS;

stoppingCriterion = STOPPING_DEFAULT;

% Initialize parameters
[m,n] = size(A) ;

tol = 1e-6 ;
tol_apg = 1e-6 ;


At = A';
G = At*A ;
opts.disp = 0;
tau = eigs(G,1,'lm',opts);
tauInv = 1/tau;

nIter = 0 ;

mu = 20 *m / norm(b,1);
muInv = 1/mu ;

lambda = ones(m,1) ;
x = zeros(n,1) ;

converged_main = 0 ;

maxIter = 200 ;
maxIter_apg = 50 ;

nz_x = (abs(x)>eps*10);

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
        case 'lambda'
            warning('The parameter LAMBDA is ignored.');
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = nan(1,maxIter) ;
errorSteps = nan(1,maxIter) ;

f = norm(x,1);
while ~converged_main
    
    lambdaScaled = muInv*lambda ;
    
    nIter = nIter + 1 ;
    
    if DISPLAY
        disp(['Iteration ' num2str(nIter) ' ||x||_0 ' num2str(sum(double(abs(x)>0)))]) ;
    end
    
    x_old_main = x ;
    
    temp = b + lambdaScaled ;
    temp = At*temp ;
    
    converged_apg = 0 ;
    
    nIter_apg = 0 ;
    
    t1 = 1 ; z = x ;
    
    muTauInv = muInv*tauInv ;
    
    while ~converged_apg
        
        nIter_apg = nIter_apg + 1 ;
        
        x_old_apg = x ;
        
        temp1 = z - tauInv*(G*z - temp) ;
        
        x = shrink(temp1, muTauInv) ;
        
        if norm(x_old_apg - x) < tol_apg*norm(x_old_apg)
            converged_apg = 1 ;
        end
        
        if nIter_apg >= maxIter_apg
            converged_apg = 1 ;
        end
        
        t2 = (1+sqrt(1+4*t1*t1))/2 ;
        z = x + ((t1-1)/t2)*(x-x_old_apg) ;
        t1 = t2 ;
        
    end
    
    lambda = lambda + mu*(b - A*x) ;
    
    timeSteps(nIter) = toc(t0) ;
    errorSteps(nIter) = norm(x-xG) ;
    
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
            nz_x = (abs(x)>eps*10);
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
            f = norm(x,1);
            criterionObjective = abs(f-prev_f)/(prev_f);
            converged_main =  ~(criterionObjective > tol);
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for ALM.');
        case STOPPING_INCREMENTS
            if norm(x_old_main - x) < tol*norm(x_old_main)
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

end

function Y = shrink(X, alpha)
    
Y = sign(X).*max(abs(X)-alpha,0) ;

end