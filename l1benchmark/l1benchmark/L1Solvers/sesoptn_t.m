function [x,diff_x, tt, ids, report] =sesoptn_t(x, x00, recData, t00, func_u, func_x, multA, multAadj, options, par)
%Sequential subspace optimization (SESOP) combined with Truncated Newton
%
%                           minimize  func_u(Ax) +func_x(x)
%
% Call:   [x,report] =sesoptn(x, func_u, func_x, multA, multAadj, options,
% par)
%
%Input:
%  x - starting point
%  func_u, func_x, multA, multAadj  - references to user functions (see below)
%  options - structure of optimization settings, provided by sesoptn_optionset.m
%
%  par  - structure of parameters, which will be transmitted to the user  functions
%  par.flagXnew can be  modified by SESOPtn:
%    1 - x has been changed from the previous call of user functions func_u and  func_x
%    0 - x the same as in previous call; persistent precomputed variables in func_u and  func_x are valid for use
%
%Output:
%   x       - optimal solution (approximation)
%   report - structure with details of optimization process
%
%User functions:
%    [f,g,HZ]=func_u(u,Z,par)
%        Input:  u=Ax
%                Z  - matrix to be multiplied by the Hessian (if needed)
%        Output: f  - function value
%                g  - gradient
%                HZ - product of the Hessian by the  matrix Z
%
%     [f,g,HW,diagH]=func_x(x,W,par);
%         Input:  x - vector argument
%                 W - matrix to be multiplied by the Hessian (if needed)
%
%     Y=multA(X,par);      % Multiply matrix X by user matrix A
%
%     Z=multAadj(Y,par);    % Multiply matrix Y by adjoint of A
%
%     If precondtioning used d_new= -Wg, then user function options.mult_precond:
%     d_new=options.mult_precond (-g, x, Ax, InsideCGforTN, par);
%            InsideCGforTN=1, when precinditioning in  process of CG for solving Newton system



% Michael Zibulevsky,
% 25.11.2007;   20.12.2007;  01.16.2008;  24.02.08; 31.03.2008; 22.06.2008; 30.06.2008; 07.07.2008; 10.07.2008;
% 13.07.2008; 29.07.2008; 05.08.2008; 08.08.2008; 03.02.2009; 16.06.2009
%
% Copyright (c) 2008 - 2009. All rights reserved. Free for academic use. No warranty

%mlock;  % Lock the function in memory to preservr persistent variables
%munlock; % unlock the function in memory to edit changes
% In order to unlock, you must stop with debugger inside this
% function and execute munlock!!

k = length(x00);

if isempty(multA),    multA     = @(x,par) x;end
if isempty(multAadj), multAadj  = @(x,par) x;end
if isempty(func_u),  func_u  = @dummy_func_u;end


global GlobalNgrads GlobalGradNorms
global GlobalNfuncs GlobalFuncValues
global GlobalNmultA
global GlobalTimeMultA GlobalTimeMultAt GlobalTimeFunc_x

% persistent sesop_iter0 u p Ap Hp P R indP indR d_tn Ad_tn normdx dx u_old x_old...
%    SubspaceGradNorms mmu ff times time0 sesop_figure_handle ...
%    timesPrecond timesA timesAt timesFunc_x plottimes timePrecond plottime ...
%    y NesWeight NemWeight d2Nem


if options.PureTN,  % 1 - Use "Classic" Truncated Newton without SESOP; 0 - TN with SESOP (recommended)
	options.nLastSteps=0;
end

if options.ShowSesopPlots,
	sesop_figure_handle=figure('Position',[800 200 600 600],'Name',options.sesop_figure_name);
end

if ~options.ContinueOldRun,

	GlobalNgrads=0; GlobalGradNorms=[];
	GlobalNfuncs=0; GlobalFuncValues=[];
	GlobalNmultA=0;
	GlobalTimeMultA=0; GlobalTimeMultAt=0; GlobalTimeFunc_x=0;
	plottime=0;
	timesPrecond=[];timePrecond=0;
	par.flagXnew=1;   % 1- compute new f,g,h in func_u and func_x (don't use old persistent f,g, Hess preparations)

	SubspaceGradNorms=[];
	mmu=[];
	ff=[];
	report.gradnorms=[];
	times=[];timesA=[];timesAt=[];timesFunc_x=[];plottimes=[];
	time0=cputime;

	sesop_iter0=0;

	p=0; %Previous step
	Ap=0;
	Hp=0;

	x0=x;

	u=multA(x,par); u_old=u;
	n=length(x);
	m=length(u);
	P=zeros(n,1);indP=0;  % Columns of P span the subspace of optimization
	R=zeros(m,1);indR=0;  % R=AP

	d_tn=0;   % Old direction for truncated newton
	Ad_tn=0;


	% f_old=1e40;  % Dummy values for stopping criterion check
	% x_old=1e40;

	if options.TESTGRAD, testgrad; return; end % Test grad, hess & adjoint operator


end

diff_x = zeros(1,options.max_sesop_iter-1);
tt = zeros(1,options.max_sesop_iter-1);
ids = zeros(1,options.max_sesop_iter-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                      Main SESOP loop
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExitFlag=0;

f_old=1e20;

for sesop_iter=sesop_iter0+1:sesop_iter0+options.max_sesop_iter,

	%cputime-time0

	if sesop_iter==1 || ~(options.FlagWolfeLineSearch && options.FlagLBFGS)
		par.flagXnew=1;
		[f_u,g_u]=func_u(u,[],par);
		[f2,g2_x]=func_x(x,[],par);
		f=f_u+f2;
		g1_x=multAadj(g_u,par);
		g_x = g1_x+g2_x;
	end

	par.flagXnew=0;
			
	%if options.FlagFISTA && sesop_iter>1
	if options.FlagFISTA && mod(sesop_iter,10)==0
		par.flagXnew=1;
		Ay=multA(y,par);
		f_u=func_u(Ay,[],par);
		f2=func_x(y,[],par);
		f=f_u+f2;
	end


	



	%%%%%%  Check stopping criteria  %%%%%%%%%%%%%%%%
    
    if sesop_iter~=1       
        diff_x(sesop_iter-1) = norm(x00 - x(1:k));
        tt(sesop_iter-1) = toc(t00);
        sci = zeros(1,recData.nSubj);
        norm_xk = norm(x(1:k),1);
        for i = 1:recData.nSubj
            sci(i) = (recData.nSubj * norm(x(((i-1)*recData.nImg+1):i*recData.nImg), 1)/norm_xk - 1) ...
                / (recData.nSubj - 1);
        end
        [dontCare, curId] = max(sci);
        ids(sesop_iter-1) = (curId == recData.id);
    end
    

        
 	if sesop_iter>sesop_iter0+1 && (norm(x-x_old)< options.dxTol || abs(f-f_old)<options.dfTol || norm(g_x)<options.gTol),
        ExitFlag=1; sesop_iter=sesop_iter-1;
        
        diff_x = diff_x(1:sesop_iter);
        tt = tt(1:sesop_iter);
        ids = ids(1:sesop_iter);
        break
    end
    
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	f_old=f;
	x_old=x;

	if   options.precond               % preconditioning
		tmp=cputime;
		d_new= -options.mult_precond (g_x, x, u, 0, par);  % New direction: precond. gradient (or PCD or SSF direction..)
		timePrecond=timePrecond+cputime-tmp;
	else
		d_new= -g_x;                                         % New direction:  gradient
	end

	%   For Debugging:
	%   [fs,gs]= fg_sesop(x,func_u,func_x,multA,multAadj,par);
	%    f, fs
	%   [g_x gs]'
	%   us=multA(x,par);
	%   u-us


	sesop_plots; %   Plot function/gradient values

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	%            Prepare SESOP directions
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if options.FlagNonlinCG

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
		%    Polak-Ribiere Nonlinear Conjugate Gradient method 
      %    (possibly with preconditioning or PCD-CG or SSF-CG)
		%    
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		gtd=g_x'*d_new;  % -d_new - [preconditioned] gradient

		if sesop_iter==1,
			P=d_new;   % P - matrix of  subspase directions for inner SESOP optimization
			%                (it is a single column now:  direction of the next CG step)
		else
			beta= (g_x'*(d_new -d_old))/gtd_old; % Polak-Ribiere CG update
			P = d_new +beta*P;
		end
		d_old=d_new;
		gtd_old=gtd;

		R=multA(P,par);

	elseif options.FlagLBFGS % MZib L-BFGS
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
		%                  Limited Memory BFGS
      %    
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
		Ldebug=0;
		if sesop_iter == 1
			P = d_new; % Initially use simple descent direction (d_new may be -g, SSF or PCD step)
			old_dgs = zeros(length(d_new),0);  % memory of dg
			old_stps = zeros(length(d_new),0); % memory of dx
			weight_Hdiag = 1;
      else
         dd=d_old-d_new;
         %weight_Hdiag=(dd'*dx)/(dd'*dd);
			[old_dgs,old_stps,weight_Hdiag] = lbfgsUpdate_mz(dd,dx,options.LBFGSmemory,Ldebug,...
				old_dgs,old_stps,weight_Hdiag);          % Update LBFGS memory
			P = lbfgs(d_new,old_dgs,old_stps,weight_Hdiag); % Compute LBFGS direction

%          dg=g_x-g_old;
% 			[old_dgs,old_stps,weight_Hdiag] = lbfgsUpdate_mz(dg,dx,options.LBFGSmemory,Ldebug,...
% 				old_dgs,old_stps,weight_Hdiag);          % Update LBFGS memory
%          Hdiag= -g_x./d_new;  % !! is valid only for smooth function!
%          P = lbfgs_mz(g_x,old_dgs,old_stps,Hdiag); % Compute LBFGS direction
%          
         
			%figure;subplot(211);imagesc(old_dgs);colorbar;subplot(212);imagesc(old_stps);colorbar
		end
		d_old=d_new;
      g_old=g_x;
      
		%if options.max_newton_iter>0, P=P/norm(P);end
      %figure;plot(P);
      if ~options.FlagWolfeLineSearch, R=multA(P,par); end

      
      %if sesop_iter>=40,keyboard;end

	elseif options.FlagNesMz % heuristic algorithm in spirit of FISTA: not in use now
		if sesop_iter==1,
			NesWeight=1;
			P= d_new;
		else
			NesWeight=0.5+sqrt(0.25+NesWeight^2);
			P= d_new +(NesWeight-1)/(NesWeight+1)*P;
		end

	elseif options.FlagFISTA 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
		% FISTA: For compatibility with sesop, we swap x and y notation 
      %        relatively to the paper of Beck and Teboulle
      %    
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		%       if sesop_iter==1,
		%          NesWeight=1
		%          y=x;
		%       else
		%          NesWeight=0.5+sqrt(0.25+NesWeight^2);
		%       end
		%       x_new=y+d_new;
		%       y=x_new+(NesWeight-1)/(NesWeight+1)*(x_new-x);
		%
		if sesop_iter==1,
			NesWeight=1
			y=x;
		else
			NesWeight=0.5+sqrt(0.25+NesWeight^2);
		end
		y_new=x+d_new;
		x=y_new+(NesWeight-1)/(NesWeight+1)*(y_new-y);
		y=y_new;

	else

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      %
		%               All other methods
      %
      %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if sesop_iter>1,  % Put previous normalized step as a new direction


			p=(1/normdx)*dx; % Direction of the last step (to be added to the subspace)

			%dx=x-(x_old+d_tn) ;       % If TN, (x_old+d_tn)  - current approx. minimizer of quadratic model

			r=(1/normdx)*(u-(u_old+Ad_tn));  %  r=A*p

			if mod(sesop_iter,options.PeriodRestoreAp)==0 %|| mod(sesop_iter+1,options.PeriodRestoreAp)==0,
				r1=multA(p,par);        % This line added  to prevent error accumulation
				%fprintf(' norm(r1-r)=%g \n',norm(r1-r));
				r=r1;
			end

			% Put a new column into matrices P and R in circular way
			if options.nLastSteps>0,
				indP=indP+1; if indP>options.nLastSteps, indP=1; end
				P(:,indP)=p; R(:,indP)=r;
			end

		end

		d_new_norm=(1/norm(d_new))*d_new;
		Ad_new_norm= multA(d_new_norm,par);

		ind_new=min(sesop_iter,options.nLastSteps+1)-1;

		if options.FlagGrad_dir,
			ind_new=ind_new+1;
			P(:,ind_new)=d_new_norm;  %Put normalized [preconditioned] gradient at the last free column
			R(:,ind_new)=Ad_new_norm;
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%   Two directions of Nemirovski: x-x0 and  sum w_k*g_k
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if options.FlagNemirovski_dir,
			if sesop_iter==1,
				NemWeight=1;
				d2Nem= -g_x;
			else
				NemWeight=0.5+sqrt(0.25+NemWeight^2);
				d2Nem=d2Nem - NemWeight*g_x;

				d1Nem=x-x0;
				Ad1Nem=multA(d1Nem,par);   % Durty way: this multiplication can be saved

				ind_new=ind_new+1;
				tmp=1/(norm(d1Nem)+1e-20);
				P(:,ind_new)=tmp*d1Nem;
				R(:,ind_new)=tmp*Ad1Nem;
			end

			Ad2Nem=multA(d2Nem,par);   % Durty way: this multiplication can be saved
			ind_new=ind_new+1;
			tmp=1/(norm(d2Nem)+1e-20);
			P(:,ind_new)=tmp*d2Nem;
			R(:,ind_new)=tmp*Ad2Nem;
		end






		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
		%
		%  Truncated Newton related directions
      %
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		if options.max_iter_CGinTN > 0,

			if sesop_iter>1,
				Ap=r;
				[f_u,g_u,h1Ap]=func_u(u,Ap,par);  % h1Ap - hess in u multiply Ap
				par.flagXnew=0;  % Only u was changed last time
				[f2,g2_x,h2p]=func_x(x,p,par);
				Hp=multAadj(h1Ap,par)+h2p;
			end


			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% In SVM  linear branches of penalty are not in use for Hessian calculation
			% Therefore we can use restricted matrix A to accelerate Hess-vector multiplication:
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			if options.FlagRestrictedMatrix_TN
				ind=find(u>-1 & u< -1+ par.eps_smooth_const_linear);
				%fprintf('Sparsity=%.4f ',length(ind)/m);
				Arestrict=par.A(ind,:);  % Data matrix restrictet to "support vectors":
				%                          those, whose Huber penalty is on quadratic branch
				multA_tn  = @(x,par)  Arestrict*x;      % TN user function  y=Ax
				multAt_tn = @(x,par)  (x'*Arestrict)';  % TN user function  y=A'*x
				func_u_tn = @func_u_restrict;
				Ad_new_norm=Ad_new_norm(ind);
				if sesop_iter>1, Ap=Ap(ind); end
			else
				multA_tn  = multA;
				multAt_tn  = multAadj;
				func_u_tn = func_u;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Solve approximately Hd=-g_x using options.max_iter_CGinTN  Conjugate Gradient iterations
			%  and taking into account the last step p:
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			[d_tn, Ad_tn, Hd_tn,  p_tn, Ap_tn, Hp_tn, d_new_tn, reportCG] = ...
				CGforSesopTN(func_u_tn, func_x, multA_tn, multAt_tn, -g_x, d_new_norm, Ad_new_norm, p, Ap,Hp, x, u, options,par);

			if 1 %options.FlagRestrictedMatrix_TN
				% tmp=zeros(m,1);tmp(ind)=Ad_tn; Ad_tn=tmp;
				% tmp=zeros(m,1);tmp(ind)=Ap_tn; Ap_tn=tmp;
				Ad_tn=multA(d_tn,par);
			end

			if ~options.PureTN      % Two directions below are not used by "Classical"  TN

				Ap_tn=multA(p_tn,par);

				ind_new=ind_new+1;
				tmp=1/(norm(p_tn)+1e-20);
				P(:,ind_new)=tmp*p_tn;
				R(:,ind_new)=tmp*Ap_tn;

				ind_new=ind_new+1;
				tmp=1/(norm(d_new_tn)+1e-20);
				P(:,ind_new)=tmp*d_new_tn;            % we normalize it mzib,09.12.2009
				R(:,ind_new)=tmp*multA(d_new_tn,par);

				ind_new=ind_new+1;
				tmp=1/(norm(d_tn)+1e-20);
				P(:,ind_new)=tmp*d_tn;
				R(:,ind_new)=tmp*Ad_tn;

			else           % Classical TN:  1d linesearch subspace P
				tmp=1/(norm(d_tn)+1e-20);
				P = tmp*d_tn;
				R = tmp*Ad_tn;
			end


		end % end Truncated Newton related directions

	end    % end  else 'if options.FlagNonlinCG'



	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	%      Newton Optimization Over Subspace or Line Search
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   u_old=u;

	if options.max_newton_iter>0,
		[x,u]=NewtonOptimizationOverSubspace(x_old,u_old,P,R);  % Nested function
		if mod(sesop_iter+1,options.PeriodRestoreAx)==0,  %sesop_iter>1 % && ~options.FlagNonlinCG
			%u1=multA(x,par);   % This line added  to prevent error accumulation
			%fprintf(' norm(u1-u)=%g \n',norm(u1-u));
			%u=u1; %fprintf('***');
			u=multA(x,par);   % This line added  to prevent error accumulation
		end

		% Debug:
		%       [f_u,g_u]=func_u(u,[],par);
		%       [f2,g2_x]=func_x(x,[],par);
		%       f=f_u+f2;
		% if f>f_old +1e-7, f, f_old, disp(''sesoptn-485: f>f_old !!!''); keyboard;end
		par.flagXnew=1;
	else  % Fixed step-size (used in SSF - separable surrogate Functions and in FISTA and Nesterov method)
	    	% or Wolfe Line Search used in LBFGS

		if options.FlagWolfeLineSearch && options.FlagLBFGS
			par.flagXnew=1;
			if 0  %Wolfe line search
				%   [t,f,g,LSfunEvals] = WolfeLineSearch(x,t,d,f,g,  gtd,   c1,  c2,LS,25,tolX, debug,doPlot,1,funObj,varargin{:});
				[t,f,g_x,LSfunEvals] =  WolfeLineSearch_mz(x,1,P,f,g_x,g_x'*P,1e-4,0.9,4,25,1e-14,1,    0,     1,@fg,     par);
			else  % Back-tracking line search
				t=1;f0=f;
				for LSfunEvals=1:40,
					[f,g_x]=fg(x+t*P,par);
					if f<f0 + 1e-6*max(abs(f0),1), break;end
					t=0.25*t;
				end
				%fprintf('%d %g;  ',LSfunEvals,t);
				x=x+t*P;
			end
			par.flagXnew=0;
      else % Fixed step size
			
			if options.FlagNesMz
				x=x+P;
			elseif ~options.FlagFISTA
				x=x+d_new;
			end
			u=multA(x,par);
			par.flagXnew=1;
		end
		
	end


	
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
	
	dx=x-(x_old+d_tn) ;       % If TN, (x_old+d_tn)  - current approx. minimizer of quadratic model
	normdx=norm(dx)+1e-30;
	
		
end %%%%%%%%%%%%% End of SESOP main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




ExitFlag=1;
sesop_plots; %   Plot function/gradient values

% if  isfield(options, 'report_func')  & ~isempty(options.report_func),
% 	options.report_func(x,report,par);                 %%%%%     User report function
% end
fprintf('\n');

report.gradnorms=GlobalGradNorms;
report.Niter=sesop_iter;
report.func_values=ff;
report.times=times;
report.timesA=timesA;
report.timesAt=timesAt;
report.timesFunc_x=timesFunc_x;
report.timesPrecond=timesPrecond;
sesop_iter0=sesop_iter;




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%               NESTED  FUNCTIONS in  SESOPTN
%                                                      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   testgrad:  Test gradient, Hessian and adjoint operator in sesoptn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	function testgrad

		% Test adjoint operator
		x=rand(size(x)); Ax=multA(x,par);
		y=rand(size(Ax)); Aadj_y=multAadj(y,par);
		tmp1=Ax(:)'*y(:); tmp2=x(:)'*Aadj_y(:);
		fprintf('<Ax,y>=%g <x,Aadj_y>=%g ratio=%g \n', tmp1,tmp2, tmp1/tmp2);

		% Test grad & Hess
		u=rand(m,1);
		[f,g_u,H_u]=func_u(u,eye(m),par);
		[g_u_num, H_u_num]=ghnum(u,func_u,[],par);
		figure; subplot(121); plot([g_u g_u_num]); title('[g_u g_u_num]');
		subplot(122); plot([g_u-g_u_num]); title('g_u-g_u_num');
		figure; subplot(121); imagesc(H_u);colorbar; title('H_u');
		subplot(122); imagesc(H_u-H_u_num);colorbar; title('H_u-H_u_num');


		x=rand(n,1);
		[f,g_x,H_x]=func_x(x,eye(n),par);
		[g_x_num, H_x_num]=ghnum(x,func_x,[],par);
		figure; subplot(121); plot([g_x g_x_num]); title('[g_x g_x_num]');
		subplot(122); plot([g_x-g_x_num]); title('g_x-g_x_num');
		figure; subplot(121); imagesc(H_x);colorbar; title('H_x');
		subplot(122); imagesc(H_x-H_x_num);colorbar; title('H_x-H_x_num');

		% 	[f,g]= fg_sesop(x,func_u,func_x,multA,multAadj,par);
		% 	[g_num]=ghnum(x, @fg_sesop,func_u,func_x,multA,multAadj,par);
		%     figure; subplot(121); plot([g g_num]); title('[SESOP total: g  g_num]');
		%     subplot(122); plot([g-g_num]); title('g - g_num');

		drawnow
	end %testgrad1
	
	
	
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       sesop_plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		
   function sesop_plots %  Plot function/gradient values
      plottime0=cputime;
      
      GlobalGradNorms=[GlobalGradNorms norm(g_x)];
      ff=[ff f]; GlobalFuncValues=ff;
      times=[times cputime-time0 - plottime];
      timesA=[timesA GlobalTimeMultA];
      timesAt=[timesAt GlobalTimeMultAt];
      timesFunc_x=[timesFunc_x GlobalTimeFunc_x];
      timesPrecond=[timesPrecond timePrecond];


      if (sesop_iter==1 || mod(sesop_iter,options.period_show_progress)==0 || ExitFlag) && options.ShowSesopPlots

         nncg=[0:length(ff)-1]*(options.max_iter_CGinTN +1);

         figure(sesop_figure_handle);
         subplot(221); semilogy(nncg,GlobalGradNorms);grid
         title('Norm of gradient with SESOP iterations');
         subplot(223); %plot(ff);
         indf=[1: max(length(ff)-1,1)];
         semilogy(nncg(indf), ff(indf) - min(ff)+1e-15); grid;
         title(sprintf('FuncValue minus Best (%g)',min(ff)));
         if sesop_iter>1 
            if options.max_newton_iter>0,
               subplot(222);semilogy(SubspaceGradNorms);grid
               title('Norm of subspace gradient with Newton iterations');
            end
            subplot(224);plot([times' timesA' timesAt' timesFunc_x' timesPrecond'...
               timesA'+timesAt'+timesFunc_x'+timesPrecond']);
            legend({'cpu time' 'time A' 'time At' 'time F(x)' 'Precond' 'sumAAtFprc'});grid
            drawnow
         end
      end



      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %      User plot of intermediate results: options.report_func()
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      report.gradnorms=GlobalGradNorms;
      report.Niter=sesop_iter-1+ExitFlag;
      report.func_values=ff;
            
      if (sesop_iter==1 || mod(sesop_iter,options.period_show_progress)==0) ...
            && isfield(options, 'report_func')  &&  ~isempty(options.report_func),
         options.report_func(x,report,par);
         par.flagXnew=1;   % 1- compute new f,g,h in func_u and func_x (don't use old persistent f,g, Hess preparations)
      end
      
      plottime=plottime+cputime-plottime0;
   end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        NewtonOptimizationOverSubspace
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





	function [x,u]=NewtonOptimizationOverSubspace(x0nwt,Ax0nwt,P,R)
		%function [x,u,SubspaceGradNorms]=NewtonOptimizationOverSubspace(x0nwt, Ax0nwt, P, R, func_u, func_x, SubspaceGradNorms, options, par)

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		%           Newton optimization min_alpha f(x+P*alpha) 
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		%  matrix R is equal to AP
		
		dim_subspace=size(P,2);
		alpha=zeros(dim_subspace,1);

		x=x0nwt;
		u=Ax0nwt;
		
		
		for i_newton=1:options.max_newton_iter

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%          Compute Newton direction
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         [f_u,g_u,hR]=func_u(u,R,par);  % hR - hess in u multiply R

			if options.AnalyticHessXmult
				[f2,g2_x,hP]=func_x(x,P,par);
			else  % Hess mult. via finite differences
				eps_hess=1e-8;
				[np,mp]=size(P);
				hP=zeros(np,mp);
				[f2,g2_x]=func_x(x,[],par);
				for ihess=1:mp
					x1=x+eps_hess*P(:,ihess);
					[f21,g2_x1]=func_x(x1,[],par);
					hP(:,ihess)=(1/eps_hess)*(g2_x1-g2_x);
				end
			end
					
					
			par.flagXnew=0;
			f=f_u+f2;
			g_alpha=(g_u'*R)' + (g2_x'*P)'; % compute g_alpha= R'*g_u + P'*g2_x in efficient way
			if i_newton==1,  SubspaceGradNorms=[SubspaceGradNorms norm(g_alpha)];  end
			h_alpha=R'*hR + P'*hP;  % Hessian in alpha

			%%% Solve modified Newton system: Invert sign of negative eigenvalues %%%%%

			[S,Lambda]=eig(h_alpha);
			lam=diag(Lambda);
			lam=abs(lam);
			lam=max(lam,1e-12*max(lam));
			d_alpha=-S*(diag(1./lam)*(S'*g_alpha)); %Newton direction in alpha
			abs_g0_alpha_d=abs(g_alpha'*d_alpha);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Back-tracking line search
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			%     ff=[];mmu=linspace(-3,3,1000);
			%     for alpha=aalpha, u=Ax0nwt+R*alpha;ff=[ff func(u,[],par)];end
			%     figure(555);plot(aalpha,ff)

			f0=f;
			step_size=1; % step size

			for j=1:40,
				alpha_new=alpha+step_size*d_alpha;
				u=Ax0nwt+R*alpha_new;
				x=x0nwt + P*alpha_new;
				par.flagXnew=1;

				[f_u,g_u]=func_u(u,[],par);
				[f2,g2_x]=func_x(x,[],par);
				par.flagXnew=0;
				f=f_u+f2;
				g_alpha=(g_u'*R)' + (g2_x'*P)';

				%if f<f0 + 1e-7*max(abs(f0),1) || g_alpha'*d_alpha<=0.1*abs_g0_alpha_d, break;end
				if f<f0 || (g_alpha'*d_alpha<= 0.1*abs_g0_alpha_d &&   f<f0 + 1e-8*max(abs(f0),1)) , break;end
				step_size=0.25*step_size;
         end
 			alpha=alpha_new;
			mmu=[mmu step_size]; %fprintf('inwt%d step_size=%g f=%g\n',i_newton,step_size,f)
			SubspaceGradNorms=[SubspaceGradNorms norm(g_alpha)+1e-20];
         %fprintf('%g ',step_size);

			if norm(g_alpha)<1e-10, break; end
      end

      %if f>f_old +1e-7, f, f_old, disp('sesoptn-690: f>f_old !!!'); keyboard;end

	end % NewtonOptimizationOverSubspace


end %%%%%%   sesopmz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                   EXTERNAL FUNCTIONS: put_column,  CGforSesopTN, ...
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















function [P, indP]=put_column(P,indP,pnew,nmax)
% Put a new column into matrix P in circular way

if nmax>0,
	indP=indP+1;
	if indP>nmax, indP=1;end
	P(:,indP)=pnew;
else
	indP=0;
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CGforSesopTN:   Solve equation Hx = b with Conjugate Gradients, starting with x=0  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [x,Ax,Hx, p,Ap,Hp,d_new,report] = ...
	CGforSesopTN(func_u, func_x, multA, multAadj, b, d_new, Ad_new, p, Ap, Hp, x_nonlin, Ax_nonlin, options,par,varargin)
%
% function [x,Hx, p,Hp,r,Hr,Ar, report] = ...
%	CGforSesopTN(func_u, func_x, multA, multAadj, b, Ab, Hb,  p, Ap, Hp, x_nonlin, Ax_nonlin, options,par,varargin)
%
%  Solve equation Hx = b with Conjugate Gradients, starting with x=0,
%  where H - symmetric, positive semidefinite matrix
%  Hx is computed using  func_u, func_x, multA, multAadj
%  p - direction of the last SESOP or CG step
%  x_nonlin - current vector of variables optimized by SESOP

%global Niter

epsilon = 1e-9; % Required relative accuracy improvement in norm of residual r=Hx-b
epsilonsq=epsilon^2;

%beta=0;
%x=x0;
%p=0;
%Hx= multH(x,par,varargin{:});
%r= Hx - b;

x=0;
Hx=0;
Ax=0;

r= -b;  % residual (also the gradient of the minimized quadratic function)

normr=norm(r(:));
sqnormr=normr^2;
sqnormr0=sqnormr;

report.gradnorms=[];
for Niter=1:options.max_iter_CGinTN,

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%          Progress report
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	report.Niter=Niter;
	report.gradnorms=[report.gradnorms normr];
	if mod(Niter,options.period_show_progressCG)==0,
		if isfield(options, 'report_funcCG'),
			options.report_funcCG(x,report,par,varargin{:});
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% CG iteration: Instead of using standard CG update formula, we minimize 
	%   explicitly quadratic function over 2D subspace spanned by previous 
	%   direction and current [preconditioned] gradient. 
	%
	%   In the case of non-positive definite 2D Hessian, we exit CG process
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%% Compute Ad_new and Hd_new    %[Hr, Ar]=multH(r, par,varargin{:});%%%%%%%
	if Niter>1, Ad_new=multA(d_new,par);end
	[f_u,g_u,h1Ad_new]=func_u(Ax_nonlin,Ad_new,par);  % h1Ad_new - hess in u multiply Ad_new
	[f2,g2_x,h2d_new]=func_x(x_nonlin,d_new,par);
	Hd_new=multAadj(h1Ad_new,par)+h2d_new;
	% par.flagXnew=0; % Next time use old f,g and Hess preparations in func_u and func_x
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if (Niter == 1) &&  ( all(p==0)  ||  options.PureTN),  % When the first iter of SESOP or "Classical" TN flag
		P=d_new;
		HP= Hd_new; 
		AP=Ad_new;
	else
		P(:,1)=p; P(:,2)= d_new; 
		HP(:,1)=Hp; HP(:,2)= Hd_new; 	
		AP(:,1)=Ap; AP(:,2)= Ad_new; 
	end
	Hsubspace= P'*HP;  % Hessian in subspace of columns of P
	g_subspace=P'*r;   % Gradient in this subspace

	%%% Solve modified Newton system: Invert sign of negative eigenvalues %%%%%

	[S,Lambda]=eig(Hsubspace);
	lam=diag(Lambda);
	if  any(lam < -1e-8*max(lam)), 
		fprintf('\n N_iterCG=%d Early exit in CGforSesopTN: negative curvature     \n',Niter);
		% Alternatively we coulld add a prox term like -1/2*min(lam)*x'*x to the objective in spirit of Trust-Region algorithm
		% More accurately, norms of p and r should be also accounted. Must consider this option in the future!
		if Niter==1,     x=d_new;  Ax=Ad_new;  Hx=Hd_new;    end
		break;
	end 
	lam=abs(lam);
	lam=max(lam,1e-12*max(lam));
	d_newton=-S*( 1 ./ lam .*(S'*g_subspace)); %Newton direction in the subspace of [p,r]
	p=P*d_newton;  
	Ap=AP*d_newton;
	Hp=HP*d_newton;
	x=x+ p;   
	Ax=Ax+Ap;
	r=r+ Hp; 
	
	if   options.precond               % preconditioning
		d_new=options.mult_precond (-r, x_nonlin, Ax_nonlin,1,par);   % New direction: precond. gradient
	else
		d_new= -r;          % New direction:  gradient (residual)
	end
	d_new=(1/norm(d_new))*d_new;

	sqnormr=r(:)'*r(:);
	if sqnormr/sqnormr0<epsilonsq, break;end

end
Hx=r+b;
if   Niter>15,
	Ax=multA(x,par); Ap=multA(p,par);  % Recompute Ax and Ap to avoid error accumulation
end

end  % function CGforSesopTN





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      func_u_restrict  (for Huber-like  penalty, used by SVM Truncated Newton)
%
%      Compute Hessian multiplication HZ only for the quadratic branch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function	[f_u,g_u,HZ]=func_u_restrict(u,Z,par)
mu= par.c;
eps=par.eps_smooth_const_linear;
HZ=(mu/eps)*Z;
f_u=[];
g_u=[];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ghnum: gradient and Hessian  approximations using forward differences 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function [gnum,hnum]=ghnum(x,funname,varargin)
%
% Compute gradient and Hessian  approximations using forward differences 
%
%  Michael Zibulevsky

eps=1e-8;
n=length(x(:));
gnum=zeros(n,1);
[f0,g0]=funname(x,varargin{:});

for i=1:n,
  x(i)=x(i)+eps;
  [f,g]=funname(x,varargin{:});
  gnum(i)=(1/eps)*(f-f0);
  if nargout >1, hnum(:,i)=(1/eps)*(g-g0); end
  x(i)=x(i)-eps;
  % x(i)=x(i)+1i*eps; % For complex argument
  % [f,g]=feval(funname,x,varargin{:});
  % gnum(i)=gnum(i)+ (1i/eps)*(f-f0);
  % x(i)=x(i)-1i*eps;
end
end % function ghnum








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              muldm  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function B=muldm(d,A)  % Compute efficiently   diag(d)*A
%B=muldm(d,A)
%
% B= diag(d)*A
d=d(:);
[M,N]=size(A);
n=length(d);
if n ~= M, error('Wrong matrix sizes');end

B=A;
for j=1:N
	B(:,j)=d.*B(:,j);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              dummy_func_u  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g,Hz]= dummy_func_u(u,Z,par);
f=0;
if nargout>1, g=0*u;end
if nargout>2, Hz=0*Z;end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

