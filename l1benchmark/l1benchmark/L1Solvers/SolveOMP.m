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

function [x, k] = SolveOMP(A, y, varargin)
% SolveOMP: Orthogonal Matching Pursuit
% Usage
%	[sols, iters, activationHist] = SolveOMP(A, y, N, maxIters, lambdaStop, solFreq, verbose, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%	maxIters    maximum number of iterations to perform. If not
%               specified, runs to stopping condition (default)
%   lambdaStop  If specified, the algorithm stops when the last coefficient 
%               entered has residual correlation <= lambdaStop.
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default)
%	OptTol      Error tolerance, default 1e-5
% Outputs
%	 sols            solution(s) of OMP
%    iters           number of iterations performed
% Description
%   SolveOMP is a greedy algorithm to estimate the solution 
%   of the sparse approximation problem
%      min ||x||_0 s.t. A*x = b
%   The implementation implicitly factors the active set matrix A(:,I)
%   using Cholesky updates. 
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as an m-file. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
% See Also
%   SolveLasso, SolveBP, SolveStOMP
%
DEBUG = 0;

STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_OBJECTIVE_VALUE;

stoppingCriterion = STOPPING_DEFAULT;

OptTol = 1e-5;
lambdaStop = 0;
maxIters = length(y);
[n,N]= size(A);

% Parameters for linsolve function
% Global variables for linsolve function
global opts opts_tr machPrec
opts.UT = true; 
opts_tr.UT = true; opts_tr.TRANSA = true;
machPrec = 1e-5;

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
            if parameterValue>maxIters
                if DEBUG>0
                    warning('Parameter maxIteration is larger than the possible value: Ignored.');
                end
            else
                maxIters = parameterValue;
            end
        case 'tolerance'
            OptTol = parameterValue;
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

% Initialize
x = zeros(N,1);
k = 1;
R_I = [];
activeSet = [];
res = y;
normy = norm(y);
resnorm = normy;
done = 0;

while ~done && k<maxIters
    corr = A'*res;
    if isNonnegative
        [maxcorr i] = max(corr);
    else
        [maxcorr i] = max(abs(corr));
    end
    
    if maxcorr<=0
        done = 1;
    else
        newIndex = i(1);
        
        % Update Cholesky factorization of A_I
        [R_I, done] = updateChol(R_I, n, N, A, activeSet, newIndex);
    end
    
    if ~done
        activeSet = [activeSet newIndex];
        
        % Solve for the least squares update: (A_I'*A_I)dx_I = corr_I
        dx = zeros(N,1);
        z = linsolve(R_I,corr(activeSet),opts_tr);
        dx(activeSet) = linsolve(R_I,z,opts);
        x(activeSet) = x(activeSet) + dx(activeSet);
        
        % Compute new residual
        res = y - A(:,activeSet) * x(activeSet);
        
        switch stoppingCriterion
            case STOPPING_SUBGRADIENT
                error('Subgradient is not a valid stopping criterion for OMP.');
            case STOPPING_DUALITY_GAP
                error('Duality gap is not a valid stopping criterion for OMP.');
            case STOPPING_SPARSE_SUPPORT
                error('Sparse support is not a valid stopping criterion for OMP.');
            case STOPPING_OBJECTIVE_VALUE
                resnorm = norm(res);
                if ((resnorm <= OptTol*normy) || ((lambdaStop > 0) && (maxcorr <= lambdaStop)))
                    done = 1;
                end
            case STOPPING_GROUND_TRUTH
                done = norm(xG-x)<OptTol;
            otherwise
                error('Undefined stopping criterion');
        end

        
        if DEBUG>0
            fprintf('Iteration %d: Adding variable %d\n', k, newIndex);
        end
        
        k = k+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, flag] = updateChol(R, n, N, A, activeSet, newIndex)
% updateChol: Updates the Cholesky factor R of the matrix 
% A(:,activeSet)'*A(:,activeSet) by adding A(:,newIndex)
% If the candidate column is in the span of the existing 
% active set, R is not updated, and flag is set to 1.

global opts_tr machPrec
flag = 0;

newVec = A(:,newIndex);

if isempty(activeSet),
    R = sqrt(sum(newVec.^2));
else
    p = linsolve(R,A(:,activeSet)'*A(:,newIndex),opts_tr);

    q = sum(newVec.^2) - sum(p.^2);
    if (q <= machPrec) % Collinear vector
        flag = 1;
    else
        R = [R p; zeros(1, size(R,2)) sqrt(q)];
    end
end

%
% Copyright (c) 2006. Yaakov Tsaig
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
