function [x, nIter, timeSteps, errorSteps] = SolveSesopPCD(A, b, varargin)

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_INCREMENTS = 5 ;
STOPPING_DEFAULT = STOPPING_DUALITY_GAP;

stoppingCriterion = STOPPING_DEFAULT;

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
        case 'lambda'
            warning('The parameter LAMBDA is ignored.');
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = [] ;
errorSteps = [] ;
t0 = tic;

%% original PCD-SESOP parameters

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 
options.max_sesop_iter  = 2000;  % Max  SESOP iterations
options.max_newton_iter=1;  % Max Newton iterations in subspace optimization (one can play with this)
options.max_iter_CGinTN=0;  % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)

options.precond=1;            % 1 - use user defined  preconditioning,  0 - don't use preconditioning

FlagPCD=1;                    %  when options.precond=1:
                              %  1 -  PCD (parallel coord. descent)
                              %  0 - diagonal precond.

options.nLastSteps=1; 
options.sesop_figure_name=sprintf('SESOPtn  %d CG steps per TN iter; FlagPCD=%d',options.max_iter_CGinTN, FlagPCD);

options.ShowSesopPlots=0;

par.weight_abs_penalty= 1e-0; % Weight of smooth abs penalty (mu)
par.eps_smooth_abs=1e-0;      % Smoothing parameter of asbsolute value approximation

par.multA= @(x,par)  multMatr(A,x);        % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);     % user function  y=A'*x

c_init = zeros(size(A,2),1);             % Starting point for optimization
par.y=b;

% Compute diag_AtA  to be  used for preconditioning
% diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),20,par); % For large matrices
diag_AtA=diag(A'*A);                                       % For small matrices


par.func_u =@ls_fgh;          % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_abs_smooth;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   User preconditioning function: d_new=mult_precond (-g, x, par)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if FlagPCD   % PCD direction 

	options.mult_precond =@(g, x, Ax, InsideCGforTN, par) mult_precond_pcd(g, x, Ax, InsideCGforTN, par,diag_AtA);

else         % Diagonal preconditioning
	
	options.mult_precond = @(g, x, Ax, InsideCGforTN, par) mult_diag_precond(g, x, Ax, InsideCGforTN, par,diag_AtA);

end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%      Perform SESOP optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nIter = 0;
for i=1:3
   par.weight_abs_penalty = 1e-2*par.weight_abs_penalty; % update weight of smooth abs penalty (mu)
   par.eps_smooth_abs     = 1e-2*par.eps_smooth_abs;      % update smoothing parameter of asbsolute value approximation
    
   [c, diff_x, temp_t]=sesoptn_t(c_init, xG, t0, par.func_u, par.func_x, par.multA, par.multAt,options,par);
  
   timeSteps = [timeSteps temp_t];
   errorSteps = [errorSteps diff_x] ;
   nIter = nIter + length(diff_x);
   if timeSteps(nIter) > maxTime
       break;
   end
   c_init=c;
   % disp(['relative error: ' num2str(norm(c00-c)/norm(c00))]);

end
x = c;