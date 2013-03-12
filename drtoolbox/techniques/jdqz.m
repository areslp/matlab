function varargout=jdqz(varargin)
%JDQZ computes a partial generalized Schur decomposition (or QZ
%  decomposition) of a pair of square matrices or operators.
%  
%  LAMBDA=JDQZ(A,B) and JDQZ(A,B) return K eigenvalues of the matrix pair
%  (A,B), where K=min(5,N) and N=size(A,1) if K has not been specified.
%  
%  [X,JORDAN]=JDQZ(A,B) returns the eigenvectors X and the Jordan
%  structure JORDAN:  A*X=B*X*JORDAN. The diagonal of JORDAN contains the
%  eigenvalues: LAMBDA=DIAG(JORDAN). JORDAN is an K by K matrix with the
%  eigenvalues on the diagonal and zero or one on the first upper diagonal
%  elements. The other entries are zero.
%  
%  [X,JORDAN,HISTORY]=JDQZ(A,B) returns also the convergence history.
%  
%  [X,JORDAN,Q,Z,S,T,HISTORY]=JDQZ(A,B) 
%  If between four and seven output arguments are required, then Q and Z
%  are N by K orthonormal, S and T are K by K upper triangular such that
%  they form a partial generalized Schur decomposition: A*Q=Z*S and
%  B*Q=Z*T. Then LAMBDA=DIAG(S)./DIAG(T) and X=Q*Y with Y the eigenvectors
%  of the pair (S,T): S*Y=T*Y*JORDAN (see also OPTIONS.Schur).
%  
%  JDQZ(A,B) 
%  JDQZ('Afun','Bfun')
%  The first input argument is either a square matrix (which can be full
%  or sparse, symmetric or nonsymmetric, real or complex), or a string
%  containing the name of an M-file which applies a linear operator to the
%  columns of a given matrix. In the latter case, the M-file, say Afun.m,
%  must return the dimension N of the problem with N = Afun([],'dimension').
%  For example, JDQZ('fft',...) is much faster than JDQZ(F,...), where F is
%  the explicit FFT matrix.
%  If another input argument is a square N by N matrix or the name of an
%  M-file, then B is this argument (regardless whether A is an M-file or a
%  matrix). If B has not been specified, then B is assumed to be the
%  identity unless A is an M-file with two output vectors of dimension N
%  with [AV,BV]=Afun(V), or with AV=Afun(V,'A') and BV=Afun(V,'B').
%  
%  The remaining input arguments are optional and can be given in
%  practically any order:
%  
%  [X,JORDAN,Q,Z,S,T,HISTORY] = JDQZ(A,B,K,SIGMA,OPTIONS)
%  [X,JORDAN,Q,Z,S,T,HISTORY] = JDQZ('Afun','Bfun',K,SIGMA,OPTIONS)
%  
%  where
%  
%      K         an integer, the number of desired eigenvalues.
%      SIGMA     a scalar shift or a two letter string.
%      OPTIONS   a structure containing additional parameters.
%  
%  If K is not specified, then K = MIN(N,5) eigenvalues are computed.
%  
%  If SIGMA is not specified, then the Kth eigenvalues largest in
%  magnitude are computed. If SIGMA is a real or complex scalar, then the
%  Kth eigenvalues nearest SIGMA are computed. If SIGMA is column vector
%  of size (L,1), then the Jth eigenvalue nearest to SIGMA(MIN(J,L))
%  is computed for J=1:K. SIGMA is the "target" for the desired eigenvalues.
%  If SIGMA is one of the following strings, then it specifies the desired 
%  eigenvalues.
%  
%    SIGMA            Specified eigenvalues
%  
%    'LM'             Largest Magnitude  
%    'SM'             Smallest Magnitude (same as SIGMA = 0)
%    'LR'             Largest Real part
%    'SR'             Smallest Real part
%    'BE'             Both Ends. Computes K/2 eigenvalues
%                     from each end of the spectrum (one more
%                     from the high end if K is odd.)
%  
%  If 'TestSpace' is 'Harmonic' (see OPTIONS), then SIGMA = 0 is the
%  default, otherwise SIGMA = 'LM' is the default.
%  
%  
%  The OPTIONS structure specifies certain parameters in the algorithm.
%  
%   Field name            Parameter                             Default
%  
%   OPTIONS.Tol           Convergence tolerance:                1e-8 
%                           norm(r) <= Tol/SQRT(K)   
%   OPTIONS.jmin          Minimum dimension search subspace V   K+5
%   OPTIONS.jmax          Maximum dimension search subspace V   jmin+5
%   OPTIONS.MaxIt         Maximum number of iterations.         100
%   OPTIONS.v0            Starting space                        ones+0.1*rand
%   OPTIONS.Schur         Gives schur decomposition             'no'
%                           If 'yes', then X and JORDAN are
%                           not computed and [Q,Z,S,T,HISTORY]
%                           is the list of output arguments.
%   OPTIONS.TestSpace     Defines the test subspace W           'Harmonic'
%                           'Standard':    W=sigma*A*V+B*V
%                           'Harmonic':    W=A*V-sigma*B*V
%                           'SearchSpace': W=V
%                            W=V is justified if B is positive
%                            definite.
%   OPTIONS.Disp          Shows size of intermediate residuals  'no'
%                           and the convergence history
%   OPTIONS.NSigma        Take as target for the second and     'no'
%                           following eigenvalues, the best  
%                           approximate eigenvalues from the 
%                           test subspace.  
%   OPTIONS.Pairs         Search for conjugated eigenpairs      'no'
%   OPTIONS.LSolver       Linear solver                         'GMRES'
%   OPTIONS.LS_Tol        Residual reduction linear solver      1,0.7,0.7^2,..
%   OPTIONS.LS_MaxIt      Maximum number it.  linear solver     5
%   OPTIONS.LS_ell        ell for BiCGstab(ell)                 4
%   OPTIONS.Precond       Preconditioner  (see below)           identity.
%   OPTIONS.Type_Precond  Way of using preconditioner           'left'
%  
%  For instance
%  
%    options=struct('Tol',1.0e-8,'LSolver','BiCGstab','LS_ell',4,'Precond',M);
%  
%  changes the convergence tolerance to 1.0e-8, takes BiCGstab as linear 
%  solver, and takes M as preconditioner (for ways of defining M, see below).
%
%
%  PRECONDITIONING. The action M-inverse of the preconditioner M (an 
%  approximation of A-lamda*B) on an N-vector V can be defined in the 
%  OPTIONS
%  
%     OPTIONS.Precond
%     OPTIONS.L_Precond     same as OPTIONS.Precond
%     OPTIONS.U_Precond
%     OPTIONS.P_Precond
%
%  If no preconditioner has been specified (or is []), then M\V=V (M is
%  the identity).
%  If Precond is an N by N matrix, say, K, then
%        M\V = K\V.
%  If Precond is an N by 2*N matrix, say, K, then
%        M\V = U\L\V, where K=[L,U], and L and U are N by N matrices.
%  If Precond is a string, say, 'Mi', then
%        if Mi(V,'L') and Mi(V,'U') return N-vectors 
%               M\V = Mi(Mi(V,'L'),'U')
%        otherwise 
%               M\V = Mi(V) or M\V=Mi(V,'preconditioner').
%  Note that Precond and A can be the same string.
%  If L_Precond and U_Precond are strings, say, 'Li' and 'Ui', 
%  respectively, then
%        M\V=Ui(Li(V)).
%  If (P_precond,) L_Precond, and U_precond are N by N matrices, say, 
%  (P,) L, and U, respectively, then
%        M\V=U\L\(P*V)      (P*M=L*U)
%
%     OPTIONS.Type_Precond
%  The preconditioner can be used as explicit left preconditioner
%  ('left', default), as explicit right preconditioner ('right') or 
%  implicitly ('impl').
%  
%
%  JDQZ without input arguments returns the options and its defaults.
%

%   Gerard Sleijpen.
%   Copyright (c) 2002
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


global Qschur Zschur Sschur Tschur ...
       Operator_MVs Precond_Solves ...
       MinvZ QastMinvZ

if nargin==0
   possibilities, return,
end

%%% Read/set parameters
[n,nselect,Sigma,kappa,SCHUR,...
   jmin,jmax,tol0,maxit,V,AV,BV,TS,DISP,PAIRS,JDV0,FIX_tol,track,NSIGMA,...
   lsolver,LSpar] = ReadOptions(varargin{1:nargin});

Qschur = zeros(n,0);    Zschur=zeros(n,0);; 
MinvZ  = zeros(n,0);    QastMinvZ=zeros(0,0); 
Sschur = []; Tschur=[]; history = []; 

%%% Return if eigenvalueproblem is trivial
if n<2
  if n==1, Qschur=1; Zschur=1; [Sschur,Tschur]=MV(1); end
  if nargout == 0, Lambda=Sschur/Tschur, else
  [varargout{1:nargout}]=output(history,SCHUR,1,Sschur/Tschur); end, 
return, end

%---------- SET PARAMETERS & STRINGS FOR OUTPUT -------------------------

if     TS==0, testspace='sigma(1)''*Av+sigma(2)''*Bv';
elseif TS==1, testspace='sigma(2)*Av-sigma(1)*Bv';
elseif TS==2, testspace='v'; 
elseif TS==3, testspace='Bv';
elseif TS==4, testspace='Av';
end

String=['\r#it=%i #MV=%3i, dim(V)=%2i, |r_%2i|=%6.1e  '];

%------------------- JDQZ -----------------------------------------------

% fprintf('Scaling with kappa=%6.4g.',kappa)

k=0; nt=0; j=size(V,2); nSigma=size(Sigma,1);
it=0; extra=0; Zero=[]; target=[]; tol=tol0/sqrt(nselect);

INITIATE=1;  JDV=0; 
rKNOWN=0; EXPAND=0; USE_OLD=0; DETECTED=0;

time=clock;
if TS ~=2
while (k<nselect & it<maxit)

   %%% Initialize target, test space and interaction matrices
   if INITIATE, % set new target
      nt=min(nt+1,nSigma); sigma = Sigma(nt,:); nlit=0; lit=0;  
      if j<2
        [V,AV,BV]=Arnoldi(V,AV,BV,sigma,jmin,nselect,tol);
        rKNOWN=0; EXPAND=0; USE_OLD=0; DETECTED=0; target=[];
        j=min(jmin,n-k);
      end
      if DETECTED & NSIGMA
         [Ur,Ul,St,Tt] = SortQZ(WAV,WBV,sigma,kappa);
         y=Ur(:,1); q=V*y; Av=AV*y; Bv=BV*y; 
         [r,z,nr,theta]=Comp_rz(RepGS(Zschur,[Av,Bv],0),kappa);
         sigma=ScaleEig(theta);
         USE_OLD=NSIGMA; rKNOWN=1; lit=10;
      end   
      NEWSHIFT= 1; 
      if DETECTED & TS<2, NEWSHIFT= ~min(target==sigma); end
      target=sigma; ttarget=sigma;
      if ischar(ttarget), ttrack=0; else, ttrack=track; end
      if NEWSHIFT 
         v=V; Av=AV; Bv=BV; W=eval(testspace);
         %%% V=RepGS(Qschur,V); [AV,BV]=MV(V); %%% more stability??
         %%% W=RepGS(Zschur,eval(testspace));  %%% dangerous if sigma~lambda
         if USE_OLD, W(:,1)=V(:,1); end, 
         W=RepGS(Zschur,W); WAV=W'*AV;  WBV=W'*BV;
      end
      INITIATE=0; DETECTED=0; JDV=0;

   end % if INITIATE

   %%% Solve the preconditioned correction equation
   if rKNOWN,
      if JDV, z=W; q=V; extra=extra+1; 
         if DISP,  fprintf('  %2i-d proj.\n',k+j-1), end 
      end
      if FIX_tol*nr>1 & ~ischar(target), theta=target; else, FIX_tol=0; end
      t=SolvePCE(theta,q,z,r,lsolver,LSpar,lit); 
      nlit=nlit+1; lit=lit+1; it=it+1;
      EXPAND=1; rKNOWN=0; JDV=0;
   end % if rKNOWN
    
   %%% Expand the subspaces and the interaction matrices
   if EXPAND
      [v,zeta]=RepGS([Qschur,V],t);
      V=[V,v]; 
      [Av,Bv]=MV(v); AV=[AV,Av]; BV=[BV,Bv]; 
      w=eval(testspace); w=RepGS([Zschur,W],w);
      WAV=[WAV,W'*Av;w'*AV]; WBV=[WBV,W'*Bv;w'*BV]; W=[W,w];
      j=j+1; EXPAND=0;

      %%% Check for stagnation
      if abs(zeta(size(zeta,1),1))/norm(zeta)<0.06, JDV=JDV0; end

   end % if EXPAND
 
   %%% Solve projected eigenproblem
   if USE_OLD
      [Ur,Ul,St,Tt]=SortQZ(WAV,WBV,ttarget,kappa,(j>=jmax)*jmin,y); 
   else
      [Ur,Ul,St,Tt]=SortQZ(WAV,WBV,ttarget,kappa,(j>=jmax)*jmin); 
   end

   %%% Compute approximate eigenpair and residual
   y=Ur(:,1); q=V*y; Av=AV*y; Bv=BV*y; 
   [r,z,nr,theta]=Comp_rz(RepGS(Zschur,[Av,Bv],0),kappa); 
   %%%=== an alternative, but less stable way of computing z =====
   % beta=Tt(1,1); alpha=St(1,1); theta=[alpha,beta];
   % r=RepGS(Zschur,beta*Av-alpha*Bv,0); nr=norm(r); z=W*Ul(:,1);
   rKNOWN=1; if nr<ttrack, ttarget=ScaleEig(theta); end

         if DISP,                                  %%% display history
            fprintf(String,it,Operator_MVs,j,nlit,nr), 
         end 
         history=[history;nr,it,Operator_MVs];    %%% save history

   %%% check convergence 
   if (nr<tol)
      %%% EXPAND Schur form
      Qschur=[Qschur,q]; Zschur=[Zschur,z];
      Sschur=[[Sschur;zeros(1,k)],Zschur'*Av]; 
      Tschur=[[Tschur;zeros(1,k)],Zschur'*Bv];  Zero=[Zero,0];
      k=k+1; 
      if ischar(target), Target(k,:)=[nt,0,0];
      else, Target(k,:)=[0,target]; end
      if DISP, ShowEig(theta,target,k); end
      if (k>=nselect), break; end;
      %%% Expand preconditioned Schur matrix MinvZ=M\Zschur
      UpdateMinvZ;
      J=[2:j]; j=j-1; Ur=Ur(:,J); Ul=Ul(:,J); 
      V=V*Ur; AV=AV*Ur; BV=BV*Ur; W=W*Ul; 
      WAV=St(J,J); WBV=Tt(J,J);
  
      rKNOWN=0; DETECTED=1;  USE_OLD=0;

      %%% check for conjugate pair
      if PAIRS & (abs(imag(theta(1)/theta(2)))>tol) 
         t=ImagVector(q); % t=conj(q); t=t-q*(q'*t);
         if norm(t)>tol, t=RepGS([Qschur,V],t,0); 
            if norm(t)>200*tol
               target=ScaleEig(conj(theta));
               EXPAND=1; DETECTED=0; 
               if DISP, fprintf('--- Checking for conjugate pair ---\n'), end
            end
         end
      end
    
      INITIATE = ( j==0 & DETECTED);

   elseif DETECTED %%% To detect whether another eigenpair is accurate enough 
      INITIATE=1; 
   end % if (nr<tol)
   
   %%% restart if dim(V)> jmax
   if j==jmax
      j=jmin; J=[1:j]; 
      Ur=Ur(:,J); Ul=Ul(:,J); 
      V=V*Ur; AV=AV*Ur; BV=BV*Ur; W=W*Ul; 
      WAV=St(J,J); WBV=Tt(J,J); 
   end % if j==jmax

end % while k
end % if TS~=2


if TS==2
Q0=Qschur; ZastQ=[];
% WAV=V'*AV; WBV=V'*BV;
while (k<nselect & it<maxit)

   %%% Initialize target, test space and interaction matrices
   if INITIATE & ( nSigma>k | NSIGMA), % set new target
      nt=min(nt+1,nSigma); sigma = Sigma(nt,:); nlit=0; lit=0;                
      if j<2
        [V,AV,BV]=Arnoldi(V,AV,BV,sigma,jmin,nselect,tol);
        rKNOWN=0; EXPAND=0; USE_OLD=0; DETECTED=0; target=[];
        j=min(jmin,n-k);; 
      end
      if DETECTED & NSIGMA
         [Ur,Ul,St,Tt]=SortQZ(WAV,WBV,sigma,kappa,1); 
         q=RepGS(Zschur,V*Ur(:,1)); [Av,Bv]=MV(q); 
         [r,z,nr,theta]=Comp_rz(RepGS(Zschur,[Av,Bv],0),kappa); 
         sigma=ScaleEig(theta); 
         USE_OLD=NSIGMA; rKNOWN=1; lit=10;
      end   
      target=sigma; ttarget=sigma;
      if ischar(ttarget), ttrack=0; else, ttrack=track; end
      if ~DETECTED
         %%% additional stabilisation. May not be needed
         %%% V=RepGS(Zschur,V); [AV,BV]=MV(V); 
         %%% end add. stab.
         WAV=V'*AV; WBV=V'*BV;
      end
      DETECTED=0; INITIATE=0; JDV=0; 
   end % if INITIATE

   %%% Solve the preconditioned correction equation
   if rKNOWN,
      if JDV, z=V; q=V; extra=extra+1; 
         if DISP,  fprintf('  %2i-d proj.\n',k+j-1), end 
      end
      if FIX_tol*nr>1 & ~ischar(target), theta=target; else, FIX_tol=0; end
      t=SolvePCE(theta,q,z,r,lsolver,LSpar,lit); 
      nlit=nlit+1; lit=lit+1; it=it+1;
      EXPAND=1; rKNOWN=0; JDV=0;
   end % if rKNOWN

   %%% expand the subspaces and the interaction matrices
   if EXPAND
      [v,zeta]=RepGS([Zschur,V],t); [Av,Bv]=MV(v); 
      WAV=[WAV,V'*Av;v'*AV,v'*Av]; WBV=[WBV,V'*Bv;v'*BV,v'*Bv];
      V=[V,v]; AV=[AV,Av]; BV=[BV,Bv]; 
      j=j+1; EXPAND=0;

      %%% Check for stagnation
      if abs(zeta(size(zeta,1),1))/norm(zeta)<0.06, JDV=JDV0; end 

   end % if EXPAND
 
   %%% compute approximate eigenpair
   if USE_OLD
      [Ur,Ul]=SortQZ(WAV,WBV,ttarget,kappa,(j>=jmax)*jmin,Ur(:,1)); 
   else
      [Ur,Ul]=SortQZ(WAV,WBV,ttarget,kappa,(j>=jmax)*jmin);
   end
     
   %%% Compute approximate eigenpair and residual
   q=V*Ur(:,1); Av=AV*Ur(:,1); Bv=BV*Ur(:,1);  
   [r,z,nr,theta]=Comp_rz(RepGS(Zschur,[Av,Bv],0),kappa); 
   rKNOWN=1; if nr<ttrack, ttarget=ScaleEig(theta); end

          if DISP,                                 %%% display history
             fprintf(String,it,Operator_MVs, j,nlit,nr), 
          end  
          history=[history;nr,it,Operator_MVs];   %%% save history
   
   %%% check convergence 
   if (nr<tol)
      %%% expand Schur form
      [q,a]=RepGS(Q0,q); a1=a(k+1,1); a=a(1:k,1);
      %%% ZastQ=Z'*Q0  
      Q0=[Q0,q]; %%% the final Qschur
      ZastQ=[ZastQ,Zschur'*q;z'*Q0]; Zschur=[Zschur,z]; Qschur=[Qschur,z]; 
      Sschur=[[Sschur;Zero],a1\(Zschur'*Av-[Sschur*a;0])];
      Tschur=[[Tschur;Zero],a1\(Zschur'*Bv-[Tschur*a;0])]; Zero=[Zero,0];
      k=k+1;
      if ischar(target), Target(k,:)=[nt,0,0];
      else, Target(k,:)=[0,target]; end
      if DISP, ShowEig(theta,target,k); end
      if (k>=nselect), break; end;  
      UpdateMinvZ;
      J=[2:j]; j=j-1; rKNOWN=0; DETECTED=1; 
      Ul=Ul(:,J); 
      V=V*Ul; AV=AV*Ul; BV=BV*Ul; 
      WAV=Ul'*WAV*Ul; WBV=Ul'*WBV*Ul; 
      Ul=eye(j); Ur=Ul;

      %%% check for conjugate pair
      if PAIRS & (abs(imag(theta(2)/theta(1)))>tol)
         t=ImagVector(q); 
         if norm(t)>tol,  
            %%% t perp Zschur, t in span(Q0,imag(q)) 
            t=t-Q0*(ZastQ\(Zschur'*t));
            if norm(t)>100*tol
               target=ScaleEig(conj(theta));
               EXPAND=1; DETECTED=0; USE_OLD=0; 
               if DISP, fprintf('--- Checking for conjugate pair ---\n'), end
            end
         end
      end
      INITIATE = ( j==0 & DETECTED);

   elseif DETECTED %%% To detect whether another eigenpair is accurate enough
      INITIATE=1;
   end % if (nr<tol)

   %%% restart if dim(V)> jmax
   if j==jmax
      j=jmin; J=[1:j];
      Ur=Ur(:,J); 
      V=V*Ur; AV=AV*Ur; BV=BV*Ur; 
      WAV=Ur'*WAV*Ur; WBV=Ur'*WBV*Ur; 
      Ur=eye(j); 
   end % if jmax

end % while k
Qschur=Q0;
end

time_needed=etime(clock,time);

if JDV0 & extra>0 & DISP
  fprintf('\n\n# j-dim. proj.: %2i\n\n',extra)
end

I=CheckSortSchur(Sigma,kappa); Target(1:length(I),:)=Target(I,:);

XKNOWN=0;

if nargout == 0 
   if ~DISP
   eigenvalues=diag(Sschur)./diag(Tschur)
   % Result(eigenvalues)
   return, end
else
  Jordan=[]; X=zeros(n,0);
  if SCHUR ~= 1
    if k>0
      [Z,D,Jor]=FindJordan(Sschur,Tschur,SCHUR); 
      DT=abs(diag(D)); DS=abs(diag(Jor));
      JT=find(DT<=tol & DS>tol); JS=find(DS<=tol & DT<=tol);
      msg=''; DT=~isempty(JT); DS=~isempty(JS); 
      if DT
          msg1='The eigenvalues'; msg2=sprintf(', %i',JT);
          msg=[msg1,msg2,' are numerically ''Inf'''];
      end, 
      if DS
          msg1='The pencil is numerically degenerated in the directions';
          msg2=sprintf(', %i',JS); 
          if DT, msg=[msg,sprintf('\n\n')]; end, msg=[msg,msg1,msg2,'.'];
      end, 
      if (DT | DS), warndlg(msg,'Unreliable directions'), end
      Jordan=Jor/D; X=Qschur*Z; XKNOWN=1;
    end 
  end
  [varargout{1:nargout}]=output(history,SCHUR,X,Jordan);

end

%-------------- display results -----------------------------------------
if DISP & size(history,1)>0
  rs=history(:,1); mrs=max(rs);
  if mrs>0, rs=rs+0.1*eps*mrs;
    subplot(2,1,1); t=history(:,2); 
    plot(t,log10(rs),'*-',t,log10(tol)+0*t,':')
    legend('log_{10} || r_{#it} ||_2')
    String=sprintf('The test subspace is computed as %s.',testspace);
    title(String)

    subplot(2,1,2); t=history(:,3);
    plot(t,log10(rs),'-*',t,log10(tol)+0*t,':')
    legend('log_{10} || r_{#MV} ||_2')
  

    String=sprintf('JDQZ with jmin=%g, jmax=%g, residual tolerance %g.',...
            jmin,jmax,tol); 
    title(String) 
    String=sprintf('Correction equation solved with %s.',lsolver);
    xlabel(String), 

    date=fix(clock);
    String=sprintf('%2i-%2i-%2i, %2i:%2i:%2i',date(3:-1:1),date(4:6));
    ax=axis; text(0.2*ax(1)+0.8*ax(2),1.2*ax(3)-0.2*ax(4),String)
    drawnow
  end

  Result(Sigma,Target,diag(Sschur),diag(Tschur),tol)

end

%------------------------ TEST ACCURACY ---------------------------------
if k>nselect & DISP
   fprintf('\n%i additional eigenpairs have been detected.\n',k-nselect)
end
if k<nselect & DISP
   fprintf('\nFailed to detect %i eigenpairs.\n',nselect-k)
end

if (k>0) & DISP
   Str='time_needed';                      texttest(Str,eval(Str))
   fprintf('\n%39s: %9i','Number of Operator actions',Operator_MVs)
   if Precond_Solves
   fprintf('\n%39s: %9i','Number of preconditioner solves',Precond_Solves)
   end
   if 1
   if SCHUR ~= 1 & XKNOWN
     % Str='norm(Sschur*Z-Tschur*Z*Jordan)'; texttest(Str,eval(Str),tol0)
     ok=1; eval('[AX,BX]=MV(X);','ok=0;')
     if ~ok, for j=1:size(X,2), [AX(:,j),BX(:,j)]=MV(X(:,j)); end, end
     Str='norm(AX*D-BX*Jor)';              texttest(Str,eval(Str),tol0)
   end
   ok=1; eval('[AQ,BQ]=MV(Qschur);','ok=0;')
   if ~ok, for j=1:size(Qschur,2), [AQ(:,j),BQ(:,j)]=MV(Qschur(:,j)); end, end
   if kappa == 1
     Str='norm(AQ-Zschur*Sschur)';         texttest(Str,eval(Str),tol0)
   else
     Str='norm(AQ-Zschur*Sschur)/kappa';   texttest(Str,eval(Str),tol0)
   end
   Str='norm(BQ-Zschur*Tschur)';           texttest(Str,eval(Str),tol0)
   I=eye(k);
   Str='norm(Qschur''*Qschur-I)';          texttest(Str,eval(Str))  
   Str='norm(Zschur''*Zschur-I)';          texttest(Str,eval(Str))
   nrmSschur=max(norm(Sschur),1.e-8);
   nrmTschur=max(norm(Tschur),1.e-8);
   Str='norm(tril(Sschur,-1))/nrmSschur';  texttest(Str,eval(Str))
   Str='norm(tril(Tschur,-1))/nrmTschur';  texttest(Str,eval(Str))
   end
   fprintf('\n==================================================\n')
end
if k==0
  disp('no eigenvalue could be detected with the required precision')
end

return
%%%======== END JDQZ ====================================================

%%%======================================================================
%%%======== PREPROCESSING ===============================================
%%%======================================================================

%%%======== ARNOLDI (for initial spaces) ================================        
function [V,AV,BV]=Arnoldi(v,Av,Bv,sigma,jmin,nselect,tol)
% Apply Arnoldi with M\(A*sigma(1)'+B*sigma(2)'), to construct an 
% initial search subspace
%

global Qschur

if ischar(sigma), sigma=[0,1]; end

[n,j]=size(v); k=size(Qschur,2); jmin=min(jmin,n-k);

if j==0 & k>0
  v=RepGS(Qschur,rand(n,1)); [Av,Bv]=MV(v); j=1; 
end

V=v; AV=Av; BV=Bv;
while j<jmin;
   v=[Av,Bv]*sigma';
   v0=SolvePrecond(v);
   if sigma(1)==0 & norm(v0-v)<tol, 
   %%%% then precond=I and target = 0: apply Arnoldi with A
      sigma=[1,0]; v0=Av;
   end
   v=RepGS([Qschur,V],v0); V=[V,v]; 
   [Av,Bv]=MV(v); AV=[AV,Av]; BV=[BV,Bv]; j=j+1; 
end % while

return
%%%======== END ARNOLDI =================================================

%%%======================================================================
%%%======== POSTPROCESSING ==============================================
%%%======================================================================

%%%======== SORT QZ DECOMPOSITION INTERACTION MATRICES ==================
function I=CheckSortSchur(Sigma,kappa)
% I=CheckSortSchur(Sigma)
%   Scales Qschur, Sschur, and Tschur such that diag(Tschur) in [0,1]
%   Reorders the Partial Schur decomposition such that the `eigenvalues'
%   (diag(S),diag(T)) appear in increasing chordal distance w.r.t. to
%   Sigma.  
%   If diag(T) is non-singular then Lambda=diag(S)./diag(T) are the
%   eigenvalues.

global Qschur Zschur Sschur Tschur

k=size(Sschur,1); if k==0, I=[]; return, end
% [AQ,BQ]=MV(Qschur);
%   Str='norm(AQ-Zschur*Sschur)';                texttest(Str,eval(Str))
%   Str='norm(BQ-Zschur*Tschur)';                texttest(Str,eval(Str))

%--- scale such that diag(Tschur) in [0,1] ----
[Tschur,D]=ScaleT(Tschur); Sschur=D\Sschur;

% kappa=max(norm(Sschur,inf)/norm(Tschur,inf),1);

s=diag(Sschur); t=diag(Tschur);
I=(1:k)'; l=size(Sigma,1);
for j=1:k 
  J0=(j:k)';
  J=SortEig(s(I(J0)),t(I(J0)),Sigma(min(j,l),:),kappa);
  I(J0)=I(J0(J));
end

if ~min((1:k)'==I)
   [Q,Z,Sschur,Tschur]=SwapQZ(eye(k),eye(k),Sschur,Tschur,I); 
   [Tschur,D2]=ScaleT(Tschur); Sschur=D2\Sschur;
   Qschur=Qschur*Q; Zschur=Zschur*(D*Z*D2);
else
   Zschur=Zschur*D;
end

return
%========================================================================
function [T,D]=ScaleT(T)
% scale such that diag(T) in [0,1] ----
   n=sign(diag(T)); n=n+(n==0); D=diag(n);
   T=D\T; IT=imag(T); RT=real(T);
   T=RT+IT.*(abs(IT)>eps*abs(RT))*sqrt(-1);
return
%%%======== COMPUTE SORTED JORDAN FORM ==================================
function [X,D,Jor]=FindJordan(S,T,SCHUR)
% [X,D,J]=FINDJORDAN(S,T)
%   For S and T k by k upper triangular matrices
%   FINDJORDAN computes the Jordan decomposition.
%   X is a k by k matrix of eigenvectors and principal vectors
%   D and J are k by matrices, D is diagonal, J is Jordan
%   such that S*X*D=T*X*J. (diag(D),diag(J)) are the eigenvalues.
%   If D is non-singular then Lambda=diag(J)./diag(D)
%   are the eigenvalues.

% coded by Gerard Sleijpen, May, 2002

k=size(S,1);
s=diag(S); t=diag(T); n=sign(t); n=n+(n==0); 
D=sqrt(conj(s).*s+conj(t).*t).*n;
S=diag(D)\S; T=diag(D)\T; D=diag(diag(T));

if k<1, 
  if k==0, X=[]; D=[]; Jor=[]; end
  if k==1, X=1; Jor=s; end
  return
end

tol=k*(norm(S,1)+norm(T,1))*eps;
[X,Jor,I]=PseudoJordan(S,T,tol);


if SCHUR == 0
for l=1:length(I)-1
  if I(l)<I(l+1)-1, 
    J=[I(l):I(l+1)-1];  
    [U,JJor]=JordanBlock(Jor(J,J),tol);
    X(:,J)=X(:,J)*U; Jor(J,J)=JJor;
  end
end
end

Jor=Jor+diag(diag(S)); Jor=Jor.*(abs(Jor)>tol);

return
%==================================================
function [X,Jor,J]=PseudoJordan(S,T,delta)
% Computes a pseudo-Jordan decomposition for the upper triangular 
% matrices S and T with ordered diagonal elements. 
% S*X*(diag(diag(T)))=T*X*(diag(diag(S))+Jor) 
% with X(:,i:j) orthonormal if its 
% columns span an invariant subspace of (S,T).

k=size(S,1); s=diag(S); t=diag(T); 

Jor=zeros(k); X=eye(k); J=1;

for i=2:k
  I=[1:i]; 
  C=t(i,1)*S(I,I)-s(i,1)*T(I,I); C(i,i)=norm(C,inf);
  if C(i,i)>0
    tol=delta*C(i,i);
    for j=i:-1:1 
      if j==1 | abs(C(j-1,j-1))>tol, break; end
    end
    e=zeros(i,1); e(i,1)=1;
    if j==i
      J=[J,i]; q=C\e; X(I,i)=q/norm(q);
    else
      q=X(I,j:i-1); 
      q=[C,T(I,I)*q;q',zeros(i-j)]\[e;zeros(i-j,1)];
      q=q/norm(q(I,1)); X(I,i)=q(I,1);
      Jor(j:i-1,i)=-q(i+1:2*i-j,1);
    end
  end
end
J=[J,k+1];

return
%==================================================
function [X,Jor,U]=JordanBlock(A,tol)
%  If A is nilpotent, then A*X=X*Jor with
%  Jor a Jordan block
%

k=size(A,1); Id=eye(k); 
U=Id; aa=A; j=k; jj=[]; J=1:k;

while j>0
  [u,s,v]=svd(aa); U(:,J)=U(:,J)*v;
  sigma=diag(s); delta=tol;
  J=find(sigma<delta); 
  if isempty(J),j=0; else, j=min(J)-1; end
  jj=[jj,j]; if j==0, break, end
  aa=v'*u*s; J=1:j; aa=aa(J,J);
end 
Jor=U'*A*U; Jor=Jor.*(abs(Jor)>tol);
l=length(jj); jj=[jj(l:-1:1),k];

l2=jj(2)-jj(1); J=jj(1)+(1:l2); 
JX=Id(:,J); X=Id;
for j=2:l
  l1=l2+1; l2=jj(j+1)-jj(j); 
  J2=l1:l2; J=jj(j)+(1:l2);
  JX=Jor*JX; D=diag(sqrt(diag(JX'*JX))); JX=JX/D;
  [Q,S,V]=svd(JX(J,:));
  JX=[JX,Id(:,J)*Q(:,J2)]; X(:,J)=JX;
end

J=[];
for i=1:l2
  for k=l:-1:1
    j=jj(k)+i; if j<=jj(k+1), J=[J,j]; end
  end
end

X=X(:,J); Jor=X\(Jor*X); X=U*X;
Jor=Jor.*(abs(Jor)>100*tol);

return
%%%======== END JORDAN FORM =============================================

%%%======== OUTPUT ======================================================
function varargout=output(history,SCHUR,X,Lambda)

global Qschur Zschur Sschur Tschur

if nargout == 1, varargout{1}=diag(Sschur)./diag(Tschur); return, end
if nargout > 2,  varargout{nargout}=history;        end
if nargout < 6 & SCHUR == 1
  if nargout >1, varargout{1}=Qschur; varargout{2}=Zschur; end
  if nargout >2, varargout{3}=Sschur; end
  if nargout >3, varargout{4}=Tschur; end
end

%-------------- compute eigenpairs --------------------------------------

if SCHUR ~= 1
   varargout{1}=X; varargout{2}=Lambda; 
   if nargout >3, varargout{3}=Qschur; varargout{4}=Zschur; end
   if nargout >4, varargout{5}=Sschur; end
   if nargout >5, varargout{6}=Tschur; end
end

return
%%%======================================================================
%%%======== UPDATE PRECONDITIONED SCHUR VECTORS =========================
%%%======================================================================
function   UpdateMinvZ

global Qschur Zschur MinvZ QastMinvZ

  [n,k]=size(Qschur);
  if k==1, MinvZ=zeros(n,0); QastMinvZ = []; end
  Minv_z=SolvePrecond(Zschur(:,k)); 
  QastMinvZ=[[QastMinvZ;Qschur(:,k)'*MinvZ],Qschur'*Minv_z];
  MinvZ=[MinvZ,Minv_z]; 

return
%%%======================================================================
%%%======== SOLVE CORRECTION EQUATION ===================================
%%%======================================================================

function [t,xtol]=SolvePCE(theta,q,z,r,lsolver,par,nit)

global Qschur Zschur

  Q=[Qschur,q]; Z=[Zschur,z];

  switch lsolver
    case 'exact'
         [t,xtol] = exact(theta,Q,Z,r);

    case {'gmres','cgstab','olsen'} 

      [MZ,QMZ]=FormPM(q,z); 
      %%% solve preconditioned system
      [t,xtol] = feval(lsolver,theta,Q,Z,MZ,QMZ,r,spar(par,nit));

  end
    
return
%------------------------------------------------------------------------
function [MZ,QMZ]=FormPM(q,z)
% compute vectors and matrices for skew projection

global Qschur MinvZ QastMinvZ

  Minv_z=SolvePrecond(z);
  QMZ=[QastMinvZ,Qschur'*Minv_z;q'*MinvZ,q'*Minv_z];
  MZ=[MinvZ,Minv_z];
    
return
%%%======================================================================
%%%======== LINEAR SOLVERS ==============================================
%%%======================================================================
function [x,xtol] = exact(theta,Q,Z,r)
% produces the exact solution if matrices are given
% Is only feasible for low dimensional matrices
% Only of interest for experimental purposes
%
global Operator_A Operator_B

n=size(r,1); 

if ischar(Operator_A)
   [MZ,QMZ]=FormPM(Q(:,end),Z(:,end)); 
   if n>200
     [x,xtol]=SolvePCE(theta,Q,Z,MZ,QMZ,r,'cgstab',[1.0e-10,500,4]); 
   else
     [x,xtol]=SolvePCE(theta,Q,Z,MZ,QMZ,r,'gmres',[1.0e-10,100]); 
   end
   return
end

k=size(Q,2);
Aug=[theta(2)*Operator_A-theta(1)*Operator_B,Z;Q',zeros(k,k)];
x=Aug\[r;zeros(k,1)]; x([n+1:n+k],:)=[]; xtol=1;

% L=eig(full(Aug)); plot(real(L),imag(L),'*'), pause
%%% [At,Bt]=MV(x); At=theta(2)*At-theta(1)*Bt; 
%%% xtol=norm(r-At+Z*(Z'*At))/norm(r); 

return

%%%===== Iterative methods ==============================================
function [r,xtol] = olsen(theta,Q,Z,MZ,M,r,par)
% returns the preconditioned residual as approximate solution
% May be sufficient in case of an excellent preconditioner

  r=SkewProj(Q,MZ,M,SolvePrecond(r)); xtol=0;

return
%------------------------------------------------------------------------
function [x,rnrm] = cgstab(theta,Q,Z,MZ,M,r,par)
% BiCGstab(ell) with preconditioning
% [x,rnrm] = cgstab(theta,Q,Z,MZ,M,r,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=r 
% where Atilde=(I-Z*Z)*(A-theta*B)*(I-Q*Q').
% using (I-MZ*(M\Q'))*inv(K) as preconditioner
%
% This function is specialized for use in JDQZ.
% integer nmv: number of matrix multiplications
% rnrm: relative residual norm
%
%  par=[tol,mxmv,ell] where 
%    integer m: max number of iteration steps
%    real tol: residual reduction
%
% rnrm: obtained residual reduction
%
% -- References: ETNA

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization --
%

global Precond_Type

tol=par(1); max_it=par(2); l=par(3); n=size(r,1);
rnrm=1; nmv=0;
 
if max_it < 2 | tol>=1, x=r; return, end
%%% 0 step of bicgstab eq. 1 step of bicgstab
%%% Then x is a multiple of b

TP=Precond_Type;
if TP==0, r=SkewProj(Q,MZ,M,SolvePrecond(r)); tr=r;
else, tr=RepGS(Z,r); end
rnrm=norm(r); snrm=rnrm; tol=tol*snrm;

sigma=1; omega=1; 
x=zeros(n,1); u=zeros(n,1);
J1=2:l+1; 
   
%%% HIST=[0,1];

if TP <2 %% explicit preconditioning
% -- Iteration loop
while (nmv < max_it)

   sigma=-omega*sigma;
   for j = 1:l,
      rho=tr'*r(:,j);  bet=rho/sigma;
      u=r-bet*u;
      u(:,j+1)=PreMV(theta,Q,MZ,M,u(:,j));
      sigma=tr'*u(:,j+1);  alp=rho/sigma;
      r=r-alp*u(:,2:j+1);
      r(:,j+1)=PreMV(theta,Q,MZ,M,r(:,j));
      x=x+alp*u(:,1);
      G(1,1)=r(:,1)'*r(:,1); rnrm=sqrt(G(1,1));
      if rnrm<tol, l=j; J1=2:l+1; r=r(:,1:l+1); break, end
   end
   nmv = nmv+2*l;

   for i=2:l+1 
     G(i,1:i)=r(:,i)'*r(:,1:i); G(1:i,i)=G(i,1:i)'; 
   end
   if TP, g=Z'*r; G=G-g'*g; end
   d=G(J1,1); gamma=G(J1,J1)\d;  
   rnrm=sqrt(real(G(1,1)-d'*gamma));   %%% compute norm in l-space
   %%% HIST=[HIST;[nmv,rnrm/snrm]];

   x=x+r(:,1:l)*gamma;
   if rnrm < tol, break, end     %%% sufficient accuracy. No need to update r,u
   omega=gamma(l,1); gamma=[1;-gamma];
   u=u*gamma; r=r*gamma; 
   if TP, g=g*gamma; r=r-Z*g; end

   % rnrm = norm(r); 
end

else %% implicit preconditioning

I=eye(2*l); v0=I(:,1:l); s0=I(:,l+1:2*l);
y0=zeros(2*l,1); V=zeros(n,2*l); 

while (nmv < max_it)

   sigma=-omega*sigma;
   y=y0; v=v0; s=s0;
   for j = 1:l,
      rho=tr'*r(:,j);  bet=rho/sigma;
      u=r-bet*u;
      if j>1,                 %%% collect the updates for x in l-space
         v(:,1:j-1)=s(:,1:j-1)-bet*v(:,1:j-1); 
      end
      [u(:,j+1),V(:,j)]=PreMV(theta,Q,MZ,M,u(:,j));
      sigma=tr'*u(:,j+1);  alp=rho/sigma;
      r=r-alp*u(:,2:j+1);
      if j>1, 
         s(:,1:j-1)=s(:,1:j-1)-alp*v(:,2:j); 
      end
      [r(:,j+1),V(:,l+j)]=PreMV(theta,Q,MZ,M,r(:,j));
      y=y+alp*v(:,1);  
      G(1,1)=r(:,1)'*r(:,1); rnrm=sqrt(G(1,1));
      if rnrm<tol, l=j; J1=2:l+1; s=s(:,1:l); break, end
   end
   nmv = nmv+2*l;

   for i=2:l+1 
     G(i,1:i)=r(:,i)'*r(:,1:i); G(1:i,i)=G(i,1:i)'; 
   end
   g=Z'*r; G=G-g'*g;         %%% but, do the orth to Z implicitly
   d=G(J1,1); gamma=G(J1,J1)\d;  
   rnrm=sqrt(real(G(1,1)-d'*gamma)); %%% compute norm in l-space
   x=x+V*(y+s*gamma);

   %%% HIST=[HIST;[nmv,rnrm/snrm]];

   if rnrm < tol, break, end  %%% sufficient accuracy. No need to update r,u
   omega=gamma(l,1); gamma=[1;-gamma];
   u=u*gamma; r=r*gamma; 
   g=g*gamma; r=r-Z*g;        %%% Do the orth to Z explicitly
                              %%% In exact arithmetic not needed, but
                              %%% appears to be more stable.

end
end

if TP==1, x=SkewProj(Q,MZ,M,SolvePrecond(x)); end
rnrm = rnrm/snrm;
%%% plot(HIST(:,1),log10(HIST(:,2)+eps),'*'), drawnow
return
%----------------------------------------------------------------------
function [v,rnrm] = gmres0(theta,Q,Z,MZ,M,v,par)
% GMRES
% [x,rnrm] = gmres(theta,Q,Z,MZ,M,v,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=b 
% where Atilde=(I-Z*Z)*(A-theta*B)*(I-Q*Q').
% using (I-MZ*(M\Q'))*inv(K) as preconditioner
%
% If used as implicit preconditioner then FGMRES.
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% rnrm: obtained residual reduction
%
% -- References: Saad & Schultz SISC 1986

% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 1998, Gerard Sleijpen

% -- Initialization

global Precond_Type

tol=par(1); max_it=par(2); n = size(v,1);
rnrm = 1; j=0;

if max_it < 2 | tol>=1, return, end 
%%% 0 step of gmres eq. 1 step of gmres
%%% Then x is a multiple of b
 
H = zeros(max_it +1,max_it); Rot=[ones(1,max_it);zeros(1,max_it)];

TP=Precond_Type;

TP=Precond_Type; 
if TP==0
  v=SkewProj(Q,MZ,M,SolvePrecond(v)); rho0 = norm(v); v = v/rho0;
else
  v=RepGS(Z,v); 
end

V = [v];
tol = tol * rnrm; 
y = [ rnrm ; zeros(max_it,1) ];

while (j < max_it) & (rnrm > tol),
  j=j+1;
  [v,w]=PreMV(theta,Q,MZ,M,v); 
  if TP 
    if TP == 2, W=[W,w]; end 
    v=RepGS(Z,v,0); 
  end
  [v,h] = RepGS(V,v); H(1:size(h,1),j) = h;
  V = [V, v]; 
  for i = 1:j-1,
    a = Rot(:,i);
    H(i:i+1,j) = [a'; -a(2) a(1)]*H(i:i+1,j);
  end
  J=[j, j+1];
  a=H(J,j);
  if a(2) ~= 0
     cs = norm(a); 
     a = a/cs; Rot(:,j) = a;
     H(J,j) = [cs; 0];
     y(J) = [a'; -a(2) a(1)]*y(J);
  end 
  rnrm = abs(y(j+1));
end

J=[1:j];  
if TP == 2
  v = W(:,J)*(H(J,J)\y(J));
else
  v = V(:,J)*(H(J,J)\y(J));
end

if TP==1, v=SkewProj(Q,MZ,M,SolvePrecond(v)); end

return
%%%======================================================================
function [v,rnrm] = gmres(theta,Q,Z,MZ,M,v,par)
% GMRES 
% [x,nmv,rnrm] = gmres(theta,Q,Z,MZ,M,v,par)
% Computes iteratively an approximation to the solution 
% of the linear system Q'*x = 0 and Atilde*x=r 
% where Atilde=(I-Z*Z)*(A-theta*B)*(I-Q*Q').
% using (I-MZ*(M\Q'))*inv(K) as preconditioner.
%
% If used as implicit preconditioner, then FGMRES.
%
% par=[tol,m] where
%  integer m: degree of the minimal residual polynomial
%  real tol: residual reduction
%
% nmv:  number of MV with Atilde
% rnrm: obtained residual reduction
%
% -- References: Saad
% Same as gmres0. However this variant uses MATLAB built-in functions
% slightly more efficient (see Sleijpen and van den Eshof).
%
%
% Gerard Sleijpen (sleijpen@math.uu.nl)
% Copyright (c) 2002, Gerard Sleijpen

global Precond_Type

% -- Initialization
tol=par(1); max_it=par(2); n = size(v,1);
j=0;

if max_it < 2 | tol>=1, rnrm=1; return, end 
%%% 0 step of gmres eq. 1 step of gmres
%%% Then x is a multiple of b
 
H = zeros(max_it +1,max_it); Gamma=1; rho=1;

TP=Precond_Type; 
if TP==0
  v=SkewProj(Q,MZ,M,SolvePrecond(v)); rho0 = norm(v); v = v/rho0;
else
  v=RepGS(Z,v); rho0=1;
end

V = zeros(n,0); W=zeros(n,0);
tol0 = 1/(tol*tol); 
%% HIST=1;
while (j < max_it) & (rho < tol0) 

  V=[V,v]; j=j+1;
  [v,w]=PreMV(theta,Q,MZ,M,v);
  if TP 
    if TP == 2, W=[W,w]; end 
    v=RepGS(Z,v,0); 
  end
  [v,h] = RepGS(V,v); 
  H(1:size(h,1),j)=h; gamma=H(j+1,j);
  
  if gamma==0, break %%% Lucky break-down
  else
    gamma= -Gamma*h(1:j)/gamma; 
    Gamma=[Gamma,gamma];
    rho=rho+gamma'*gamma;
  end     
          
  %% HIST=[HIST;(gamma~=0)/sqrt(rho)]; 
    
end

if gamma==0; %%% Lucky break-down
   e1=zeros(j,1); e1(1)=rho0; rnrm=0; 
   if TP == 2
     v=W*(H(1:j,1:j)\e1); 
   else
     v=V*(H(1:j,1:j)\e1); 
   end 
else %%% solve in least square sense 
   e1=zeros(j+1,1); e1(1)=rho0; rnrm=1/sqrt(rho);
   if TP == 2
     v=W*(H(1:j+1,1:j)\e1); 
   else
     v=V*(H(1:j+1,1:j)\e1); 
   end 
end

if TP==1, v=SkewProj(Q,MZ,M,SolvePrecond(v)); end
%% HIST=log10(HIST+eps); J=[0:size(HIST,1)-1]';
%% plot(J,HIST(:,1),'*'); drawnow
return
%%%======== END SOLVE CORRECTION EQUATION ===============================       

%%%======================================================================
%%%======== BASIC OPERATIONS ============================================
%%%======================================================================
function [Av,Bv]=MV(v)
% [y,z]=MV(x)
%  y=A*x, z=B*x

%  y=MV(x,theta)
%  y=(A-theta*B)*x
%

global Operator_Form Operator_MVs Operator_A Operator_B Operator_Params

  Bv=v;
  switch Operator_Form
     case 1 % both Operator_A and B are strings 
        Operator_Params{1}=v;
        Av=feval(Operator_A,Operator_Params{:}); 
        Bv=feval(Operator_B,Operator_Params{:});
     case 2
        Operator_Params{1}=v;
        [Av,Bv]=feval(Operator_A,Operator_Params{:});
     case 3
        Operator_Params{1}=v;
        Operator_Params{2}='A';
        Av=feval(Operator_A,Operator_Params{:});
        Operator_Params{2}='B';
        Bv=feval(Operator_A,Operator_Params{:});
     case 4
        Operator_Params{1}=v;
        Av=feval(Operator_A,Operator_Params{:}); 
        Bv=Operator_B*v;
     case 5
        Operator_Params{1}=v;
        Av=feval(Operator_A,Operator_Params{:});
     case 6
        Av=Operator_A*v; 
        Operator_Params{1}=v;
        Bv=feval(Operator_B,Operator_Params{:});
     case 7
        Av=Operator_A*v; 
        Bv=Operator_B*v;
     case 8
        Av=Operator_A*v;
  end

  Operator_MVs = Operator_MVs +size(v,2);

% [Av(1:5,1),Bv(1:5,1)], pause
return
%------------------------------------------------------------------------
function y=SolvePrecond(y);

global Precond_Form Precond_L Precond_U Precond_P Precond_Params Precond_Solves

    switch Precond_Form
      case 0,     
      case 1,     Precond_Params{1}=y; 
                  y=feval(Precond_L,Precond_Params{:}); 
      case 2,     Precond_Params{1}=y; Precond_Params{2}='preconditioner';
                  y=feval(Precond_L,Precond_Params{:});
      case 3,     Precond_Params{1}=y; 
                  Precond_Params{1}=feval(Precond_L,Precond_Params{:}); 
                  y=feval(Precond_U,Precond_Params{:});
      case 4,     Precond_Params{1}=y; Precond_Params{2}='L'; 
                  Precond_Params{1}=feval(Precond_L,Precond_Params{:}); 
                  Precond_Params{2}='U';
                  y=feval(Precond_L,Precond_Params{:});
      case 5,     y=Precond_L\y;
      case 6,     y=Precond_U\(Precond_L\y);
      case 7,     y=Precond_U\(Precond_L\(Precond_P*y));
    end

    if Precond_Form
      Precond_Solves = Precond_Solves +size(y,2);
    end
%% y(1:5,1), pause
return
%------------------------------------------------------------------------
function [v,u]=PreMV(theta,Q,Z,M,v)
% v=Atilde*v

global Precond_Type

  if Precond_Type
    u=SkewProj(Q,Z,M,SolvePrecond(v)); 
    [v,w]=MV(u); v=theta(2)*v-theta(1)*w;
  else
    [v,u]=MV(v); u=theta(2)*v-theta(1)*u;
    v=SkewProj(Q,Z,M,SolvePrecond(u));
  end
  
return
%------------------------------------------------------------------------
function  r=SkewProj(Q,Z,M,r);

   if ~isempty(Q), 
      r=r-Z*(M\(Q'*r));
   end 

return
%------------------------------------------------------------------------
function ppar=spar(par,nit)

k=size(par,2)-2;
ppar=par(1,k:k+2);

if k>1
   if nit>k
      ppar(1,1)=par(1,k)*((par(1,k)/par(1,k-1))^(nit-k));
   else
      ppar(1,1)=par(1,max(nit,1));
   end
end

ppar(1,1)=max(ppar(1,1),1.0e-8);

return
%------------------------------------------------------------------------
function u=ImagVector(u)
% computes "essential" imaginary part of a vector
  maxu=max(u); maxu=maxu/abs(maxu); u=imag(u/maxu);
return
%------------------------------------------------------------------------
function Sigma=ScaleEig(Sigma)
% 
%
  n=sign(Sigma(:,2)); n=n+(n==0);
  d=sqrt((Sigma.*conj(Sigma))*[1;1]).*n;
  Sigma=diag(d)\Sigma;

return
%%%======== COMPUTE r AND z =============================================
function [r,z,nrm,theta]=Comp_rz(E,kappa)
%
% [r,z,nrm,theta]=Comp_rz(E)
%   computes the direction r of the minimal residual,
%   the left projection vector z,
%   the approximate eigenvalue theta
%
% [r,z,nrm,theta]=Comp_rz(E,kappa)
%   kappa is a scaling factor.
%

% coded by Gerard Sleijpen, version Januari 7, 1998

  if nargin == 1
    kappa=norm(E(:,1))/norm(E(:,2)); kappa=2^(round(log2(kappa)));
  end

  if kappa ~=1, E(:,1)=E(:,1)/kappa; end

  [Q,sigma,u]=svd(E,0);   %%% E*u=Q*sigma, sigma(1,1)>sigma(2,2)
  r=Q(:,2); z=Q(:,1); nrm=sigma(2,2); 
  % nrm=nrm/sigma(1,1); nrmz=sigma(1,1)

  u(1,:)=u(1,:)/kappa; theta=[-u(2,2),u(1,2)];
  
return
%%%======== END computation r and z =====================================
%%%======================================================================
%%%======== Orthogonalisation ===========================================
%%%======================================================================
function [V,R]=RepGS(Z,V,gamma)
%
% Orthonormalisation using repeated Gram-Schmidt
%   with the Daniel-Gragg-Kaufman-Stewart (DGKS) criterion
%
% Q=RepGS(V)
%   The n by k matrix V is orthonormalized, that is,
%   Q is an n by k orthonormal matrix and 
%   the columns of Q span the same space as the columns of V
%   (in fact the first j columns of Q span the same space
%   as the first j columns of V for all j <= k).
%
% Q=RepGS(Z,V)
%  Assuming Z is n by l orthonormal, V is orthonormalized against Z:
%  [Z,Q]=RepGS([Z,V])
%
% Q=RepGS(Z,V,gamma)
%  With gamma=0, V is only orthogonalized against Z
%  Default gamma=1 (the same as Q=RepGS(Z,V))
%
% [Q,R]=RepGS(Z,V,gamma)
%  if gamma == 1, V=[Z,Q]*R; else, V=Z*R+Q; end
 
% coded by Gerard Sleijpen, March, 2002

% if nargin == 1, V=Z; Z=zeros(size(V,1),0); end   
if nargin <3, gamma=1; end

[n,dv]=size(V); [m,dz]=size(Z);

if gamma, l0=min(dv+dz,n); else, l0=dz; end
R=zeros(l0,dv); 

if dv==0, return, end
if dz==0 & gamma==0, return, end

% if m~=n
%   if m<n, Z=[Z;zeros(n-m,dz)]; end
%   if m>n, V=[V;zeros(m-n,dv)]; n=m; end
% end

if (dz==0 & gamma)
   j=1; l=1; J=1;
   q=V(:,1); nr=norm(q); R(1,1)=nr;
   while nr==0, q=rand(n,1); nr=norm(q); end, V(:,1)=q/nr;
   if dv==1, return, end    
else
   j=0; l=0; J=[];
end

while j<dv,
   j=j+1; q=V(:,j); nr_o=norm(q); nr=eps*nr_o;
   if dz>0, yz=Z'*q;     q=q-Z*yz;      end
   if l>0,  y=V(:,J)'*q; q=q-V(:,J)*y;  end
   nr_n=norm(q);
  
   while (nr_n<0.5*nr_o & nr_n > nr)
      if dz>0, sz=Z'*q;     q=q-Z*sz;     yz=yz+sz; end
      if l>0,  s=V(:,J)'*q; q=q-V(:,J)*s; y=y+s;    end
      nr_o=nr_n; nr_n=norm(q);                
   end
   if dz>0, R(1:dz,j)=yz; end
   if l>0,  R(dz+J,j)=y;  end

   if ~gamma
     V(:,j)=q;
   elseif l+dz<n, l=l+1; 
     if nr_n <= nr % expand with a random vector
       % if nr_n==0
       V(:,l)=RepGS([Z,V(:,J)],rand(n,1));
       % else % which can be numerical noice
       %   V(:,l)=q/nr_n;
       % end
     else
       V(:,l)=q/nr_n; R(dz+l,j)=nr_n; 
     end
     J=[1:l];
   end 

end % while j

if gamma & l<dv, V=V(:,J); end

return
%%%======== END  Orthogonalisation ======================================

%%%======================================================================
%%%======== Sorts Schur form ============================================
%%%======================================================================
function [Q,Z,S,T]=SortQZ(A,B,tau,kappa,gamma,u)
%
% [Q,Z,S,T]=SortQZ(A,B,tau)
%   A and B are k by k matrices, tau is a complex pair [alpha,beta].
%   SortQZ computes the qz-decomposition of (A,B) with prescribed
%   ordering: A*Q=Z*S, B*Q=Z*T; 
%             Q and Z are unitary k by k matrices,
%             S and T are upper triangular k by k matrices.
%   The ordering is as follows:
%   (diag(S),diag(T)) are the eigenpairs of (A,B) ordered
%   with increasing "chordal distance" w.r.t. tau.
%
%   If tau is a scalar then [tau,1] is used.
%   Default value for tau is tau=0, i.e., tau=[0,1].
%
% [Q,Z,S,T]=SortQZ(A,B,tau,kappa)
%   kappa scales A first: A/kappa. Default kappa=1.
%
% [Q,Z,S,T]=SortQZ(A,B,tau,kappa,gamma)
%   Sorts the first MAX(gamma,1) elements. Default: gamma=k
%
% [Q,Z,S,T]=SortQZ(A,B,tau,kappa,gamma,u)
%   Now, with ordering such that angle u and Q(:,1) is less than 45o,
%   and, except for the first pair, (diag(S),diag(T)) are 
%   with increasing "chordal distance" w.r.t. tau.
%   If such an ordering does not exist, ordering is as without u.
%

% coded by Gerard Sleijpen, version April, 2002

  k=size(A,1); 
  if k==1, Q=1; Z=1; S=A;T=B; return, end
  kk=k-1;
  if nargin < 3, tau=[0,1]; end 
  if nargin < 4, kappa=1;   end 
  if nargin > 4, kk=max(1,min(gamma,k-1)); end

  %%    kappa=max(norm(A,inf)/max(norm(B,inf),1.e-12),1);
  %%    kappa=2^(round(log2(kappa)));

%%%------ compute the qz factorization -------
  [S,T,Z,Q]=qz(A,B); Z=Z';
                 % kkappa=max(norm(A,inf),norm(B,inf));
                 % Str='norm(A*Q-Z*S)';texttest(Str,eval(Str))
                 % Str='norm(B*Q-Z*T)';texttest(Str,eval(Str))

%%%------ scale the eigenvalues --------------
  t=diag(T); n=sign(t); n=n+(n==0); D=diag(n);  
  Q=Q/D; S=S/D; T=T/D; 

%%%------ sort the eigenvalues ---------------
  I=SortEig(diag(S),real(diag(T)),tau,kappa);
               
%%%------ swap the qz form -------------------
  [Q,Z,S,T]=SwapQZ(Q,Z,S,T,I(1:kk)); 
                % Str='norm(A*Q-Z*S)';texttest(Str,eval(Str))
                % Str='norm(B*Q-Z*T)';texttest(Str,eval(Str))



  if nargin < 6 | size(u,2) ~= 1
     return
  else
%%% repeat SwapQZ if angle is too small
     kk=min(size(u,1),k); J=1:kk; u=u(J,1)'*Q(J,:); 
     if abs(u(1,1))>0.7, return, end
     for j=2:kk
        J=1:j;
        if norm(u(1,J))>0.7
          J0=[j,1:j-1];  
          [Qq,Zz,Ss,Tt]=SwapQZ(eye(j),eye(j),S(J,J),T(J,J),J0);
          if abs(u(1,J)*Qq(:,1))>0.71
             Q(:,J)=Q(:,J)*Qq; Z(:,J)=Z(:,J)*Zz; 
             S(J,J)=Ss; S(J,j+1:k)=Zz'*S(J,j+1:k);
             T(J,J)=Tt; T(J,j+1:k)=Zz'*T(J,j+1:k);
                 % Str='norm(A*Q-Z*S)';texttest(Str,eval(Str))
                 % Str='norm(B*Q-Z*T)';texttest(Str,eval(Str))
                 fprintf('  Took %2i:%6.4g\n',j,S(1,1)./T(1,1))
             return
          end
        end
     end
     disp(['  Selection problem:  took ',num2str(1)])
     return
  end

%%%======================================================================
function I=SortEig(s,t,sigma,kappa);
%
% I=SortEig(S,T,SIGMA) sorts the indices of [S,T] as prescribed by SIGMA
%   S and T are K-vectors.
%
%   If SIGMA=[ALPHA,BETA] is a complex pair then 
%     if CHORDALDISTANCE
%       sort [S,T] with increasing chordal distance w.r.t. SIGMA.
%     else 
%       sort S./T with increasing distance w.r.t. SIGMA(1)/SIGMA(2)
%
%  The chordal distance D between a pair A and a pair B is 
%  defined as follows.
%    Scale A by a scalar F such that NORM(F*A)=1.
%    Scale B by a scalar G such that NORM(G*B)=1.
%    Then D(A,B)=SQRT(1-ABS((F*A)*(G*B)')).
%
% I=SortEig(S,T,SIGMA,KAPPA). Kappa is a caling that effects the
% chordal distance: [S,KAPPA*T] w.r.t. [SIGMA(1),KAPPA*SIGMA(2)].

% coded by Gerard Sleijpen, version April 2002

global CHORDALDISTANCE

if ischar(sigma)
  warning off, s=s./t; warning on
  switch sigma
    case 'LM'
      [s,I]=sort(-abs(s));
    case 'SM'
      [s,I]=sort(abs(s));
    case 'LR';
      [s,I]=sort(-real(s));
    case 'SR';
      [s,I]=sort(real(s));
    case 'BE';
      [s,I]=sort(real(s)); I=twistdim(I,1);
  end
elseif CHORDALDISTANCE
  if kappa~=1, t=kappa*t; sigma(2)=kappa*sigma(2); end
  n=sqrt(s.*conj(s)+t.*t);
  [s,I]=sort(-abs([s,t]*sigma')./n);
else
  warning off, s=s./t; warning on
  if sigma(2)==0; [s,I]=sort(-abs(s));
  else, [s,I]=sort(abs(s-sigma(1)/sigma(2)));
  end 
end

return
%------------------------------------------------------------------------
function t=twistdim(t,k)

  d=size(t,k); J=1:d; J0=zeros(1,2*d);
  J0(1,2*J)=J; J0(1,2*J-1)=flipdim(J,2); I=J0(1,J);
  if k==1, t=t(I,:); else, t=t(:,I); end

return
%%%======================================================================
function [Q,Z,S,T]=SwapQZ(Q,Z,S,T,I)
% [Q,Z,S,T]=SwapQZ(QQ,ZZ,SS,TT,P)
%    QQ and ZZ are K by K unitary,  SS and TT are K by K uper triangular.
%    P is the first part of a permutation of (1:K)'.
%
%    Then Q and Z are K by K unitary, S and T are K by K upper triangular,
%    such that, for A = ZZ*SS*QQ' and B = ZZ*T*QQ', we have 
%    A*Q = Z*S, B*Q = Z*T  and LAMBDA(1:LENGTH(P))=LLAMBDA(P) where 
%    LAMBDA=DIAG(S)./DIAGg(T) and LLAMBDA=DIAG(SS)./DIAG(TT).
%
%    Computation uses Givens rotations. 
%

% coded by Gerard Sleijpen, version October 12, 1998
  

  kk=min(length(I),size(S,1)-1);
  j=1; while (j<=kk & j==I(j)), j=j+1; end
  while j<=kk
    i=I(j);
    for k = i-1:-1:j, 
      %%% i>j, move ith eigenvalue to position j 
      J = [k,k+1]; 
      q = T(k+1,k+1)*S(k,J) - S(k+1,k+1)*T(k,J);
      if q(1) ~= 0 
        q = q/norm(q);
        G = [[q(2);-q(1)],q'];
        Q(:,J) = Q(:,J)*G; 
        S(:,J) = S(:,J)*G; T(:,J) = T(:,J)*G;
      end 
      if abs(S(k+1,k))<abs(T(k+1,k)), q=T(J,k); else q=S(J,k); end
      if q(2) ~= 0
        q=q/norm(q);
        G = [q';q(2),-q(1)];
        Z(:,J) = Z(:,J)*G'; 
        S(J,:) = G*S(J,:); T(J,:) = G*T(J,:);
      end 
      T(k+1,k) = 0;
      S(k+1,k) = 0; 
    end
    I=I+(I<i); 
    j=j+1; while (j<=kk & j==I(j)), j=j+1; end
  end

return
%------------------------------------------------------------------------
function [Q,Z,S,T]=SwapQZ0(Q,Z,S,T,I)
%
% [Q,Z,S,T]=sortqz0(A,B,s,t,k)
%    A and B are k by k matrices, t and s are k vectors such that
%    (t(i),s(i)) eigenpair (A,B), i.e. t(i)*A-s(i)*B singular.
%    Computes the Schur form with a ordering prescribed by (t,s):
%       A*Q=Z*S, B*Q=Z*T such that diag(S)./diag(T)=s./t.
%    Computation uses Householder reflections. 
%

% coded by Gerard Sleijpen, version October 12, 1997
  
  k=size(S,1); s=diag(S); t=diag(T); s=s(I,1); t=t(I,1);

  for i=1:k-1
    %%% compute q s.t. C*q=(t(i,1)*S-s(i,1)*T)*q=0
    C=t(i,1)*S(i:k,i:k)-s(i,1)*T(i:k,i:k);
    [q,r,p]=qr(C);          %% C*P=Q*R
    %% check whether last but one diag. elt r nonzero
    j=k-i; while abs(r(j,j))<eps*norm(r); j=j-1; end; j=j+1;
    r(j,j)=1; e=zeros(j,1); e(j,1)=1;
    q=p*([r(1:j,1:j)\e;zeros(k-i+1-j,1)]); q=q/norm(q);%% C*q
    %%% end computation q
    z=conj(s(i,1))*S(i:k,i:k)*q+conj(t(i,1))*T(i:k,i:k)*q; z=z/norm(z);
    a=q(1,1); if a ~=0, a=abs(a)/a; q=a*q; end
    a=z(1,1); if a ~=0, a=abs(a)/a; z=a*z; end
    q(1,1)=q(1,1)+1; q=q/norm(q); q=[zeros(i-1,1);q]; 
    z(1,1)=z(1,1)+1; z=z/norm(z); z=[zeros(i-1,1);z];
    S=S-(S*q)*(2*q)'; S=S-(2*z)*(z'*S);
    T=T-(T*q)*(2*q)'; T=T-(2*z)*(z'*T);
    Q=Q-(Q*q)*(2*q)'; Z=Z-(Z*z)*(2*z)';
  end
return

%%%======== END sort QZ decomposition interaction matrices ==============

%%%======================================================================
%%%======== INITIALIZATION ==============================================
%%%======================================================================
function MyClear

global Operator_Form Operator_A Operator_B Operator_Params ...
       Precond_L Precond_U Precond_P Precond_Params ...
       Precond_Form Precond_Type ...
       Operator_MVs Precond_Solves ...
       CHORDALDISTANCE ...
       Qschur Zschur Sschur Tschur ...
       MinvZ QastMinvZ

return
%%%======================================================================
function [n,nselect,Sigma,kappa,SCHUR,...
          jmin,jmax,tol,maxit,V,AV,BV,TS,DISP,PAIRS,JDV,FIX,track,NSIGMA,...
          lsolver,par] = ReadOptions(varargin)
% Read options and set defaults

global Operator_Form Operator_A Operator_B Operator_Params ...
       Precond_Form Precond_L Precond_U Precond_P Precond_Params ...
       CHORDALDISTANCE


Operator_A = varargin{1};

n=CheckMatrix(Operator_A,1);

% defaults              %%%% search for 'xx' in fieldnames
nselect0= 5; 
maxit   = 200;          %%%% 'ma'
SCHUR   = 0;            %%%% 'sch'
tol     = 1e-8;         %%%% 'to'
DISP    = 0;            %%%% 'di'
p0      = 5; %%% jmin=nselect+p0 %%%% 'jmi'
p1      = 5; %%% jmax=jmin+p1    %%%% 'jma'
TS      = 1;            %%%% 'te'
PAIRS   = 0;            %%%% 'pai'
JDV     = 0;            %%%% 'av'
track   = 1e-4;         %%%% 'tr'
FIX     = 1000;         %%%% 'fix'
NSIGMA  = 0;            %%%% 'ns'
CHORD   = 1;            %%%% 'ch'
lsolver = 'gmres';      %%%% 'lso'
ls_maxit= 200;          %%%% 'ls_m'
ls_tol  = [0.7,0.49];   %%%% 'ls_t' 
ell     = 4;            %%%% 'ls_e'
TP      = 0;            %%%% 'ty'
L       = [];           %%%% 'l_'
U       = [];           %%%% 'u_'
P       = [];           %%%% 'p_'
kappa   = 1;            %%%% 'sca'
V0      = 'ones(n,1)+rand(n,1)';  %%%% 'v0'

%% initiation
nselect=[]; Sigma=[]; options=[]; Operator_B=[];
jmin=-1; jmax=-1; V=[]; AV=[]; BV=[]; par=[];

%------------------------------------------------
%------- Find quantities ------------------------
%------------------------------------------------
jj=[];
for j = 2:nargin
   if isstruct(varargin{j})
      options = varargin{j};
   elseif ischar(varargin{j})
      s=varargin{j}; 
      if exist(s)==2 & isempty(Operator_B)
         Operator_B=s;
      elseif  length(s) == 2 & isempty(Sigma)
        s=upper(s);
        switch s
          case {'LM','SM','LR','SR','BE'}, Sigma=s;
          otherwise
            jj=[jj,j];
        end
      else
        jj=[jj,j];
      end
   elseif min([n,n]==size(varargin{j})) & isempty(Operator_B)
      Operator_B=varargin{j}; 
   elseif length(varargin{j}) == 1
      s = varargin{j};
      if isempty(nselect) & isreal(s) & (s == fix(s)) & (s > 0)
         nselect = min(n,s);
      elseif isempty(Sigma)
         Sigma = s; 
      else
        jj=[jj,j];
      end
   elseif min(size(varargin{j}))==1  & isempty(Sigma)
      Sigma = varargin{j}; if size(Sigma,1)==1, Sigma=Sigma'; end 
   elseif min(size(varargin{j}))==2  & isempty(Sigma)
      Sigma = varargin{j}; if size(Sigma,2)>2 , Sigma=Sigma'; end 
   else
      jj=[jj,j];
   end
end

%------- find parameters for operators -----------
Operator_Params=[]; Operator_Params{2}='';
k=length(jj);
if k>0
  Operator_Params(3:k+2)=varargin(jj);
  if ~ischar(Operator_A)
    msg=sprintf(', %i',jj);
    msg=sprintf('Input argument, number%s, not recognized.',msg);
    button=questdlg(msg,'Input arguments','Ignore','Stop','Ignore');
    if strcmp(button,'Stop'), n=-1; return, end
  end
end
%------- operator B -----------------------------
if isempty(Operator_B)
   if ischar(Operator_A)
      Operator_Form=2; % or Operator_Form=3, or Operator_Form=5;
   else
      Operator_Form=8;
   end
else
   if ischar(Operator_B)
      if ischar(Operator_A), Operator_Form=1; else, Operator_Form=6; end
   elseif ischar(Operator_A)
      Operator_Form=4;
   else
      Operator_Form=7;
   end
end
if n<2, return, end

%------- number of eigs to be computed ----------
if isempty(nselect), nselect=min(n,nselect0); end

%------------------------------------------------
%------- Analyse Options ------------------------
%------------------------------------------------
fopts = [];
if ~isempty(options), fopts = fieldnames(options); end

%------- preconditioner -------------------------
Precond_L=findfield(options,fopts,'pr',[]);
[L,ok]=findfield(options,fopts,'l_',Precond_L);
if ok & ~isempty(Precond_L),
   msg =sprintf('A preconditioner is defined in');
   msg =[msg,sprintf('\n''Precond'', but also in ''L_precond''.')];
   msg=[msg,sprintf('\nWhat is the correct one?')];
   button=questdlg(msg,'Preconditioner','L_Precond','Precond','L_Precond');
   if strcmp(button,'L_Precond'), 
     Precond_L = L;
   end
else
   Precond_L = L;
end

if ~isempty(Precond_L)
  Precond_U=findfield(options,fopts,'u_',[]);
  Precond_P=findfield(options,fopts,'p_',[]);
end

Precond_Params=[]; Precond_Params{2}='';
Params=findfield(options,fopts,'par',[]);
[l,k]=size(Params);
if k>0, 
  if iscell(Params), Precond_Params(3:k+2)=Params; 
  else, Precond_Params{3}=Params; end
end
TP=findfield(options,fopts,'ty',TP);
n=SetPrecond(n,TP); if n<2, return, end

%------- max, min dimension search subspace ------
jmin=min(n,findfield(options,fopts,'jmi',jmin));
jmax=min(n,findfield(options,fopts,'jma',jmax));
if jmax < 0
   if jmin<0, jmin=min(n,nselect+p0); end
   jmax=min(n,jmin+p1); 
else
   if jmin<0, jmin=max(1,jmax-p1); end
end 
maxit=findfield(options,fopts,'ma',maxit);

%------- initial search subspace ----------------
V=findfield(options,fopts,'v',[]);
[m,d]=size(V); 
if m~=n 
  if m>n, V = V(1:n,:); end
  if m<n, V = [V;0.001*rand(n-m-1,d)]; end
end
V=orth(V); [m,d]=size(V);
if d==0, nr=0; while nr==0, V=eval(V0); nr=norm(V); V=V/nr; end, end


%------- Check definition B, Compute AV, BV -----
[AV,BV,n]=CheckDimMV(V); if n<2, return, end

%------- Other options --------------------------
tol=findfield(options,fopts,'to',tol);
kappa  = findfield(options,fopts,'sca',kappa);
kappa  = abs(kappa(1,1)); if kappa==0, kappa=1; end
PAIRS  = findfield(options,fopts,'pai',PAIRS,[0,1]);
SCHUR  = findfield(options,fopts,'sch',SCHUR,[0,1,0.5]);
DISP   = findfield(options,fopts,'di',DISP,[0,1]);
JDV    = findfield(options,fopts,'av',JDV,[0,1]);
track  = max(abs(findfield(options,fopts,'tr',track,[0,track,inf])),0);
NSIGMA = findfield(options,fopts,'ns',NSIGMA,[0,1]);
FIX    = max(abs(findfield(options,fopts,'fix',0,[0,FIX,inf])),0);
CHORDALDISTANCE = findfield(options,fopts,'ch',CHORD,[0,1]);

[TS0,ok] = findfield(options,fopts,'te',TS);

if ok & ischar(TS0)
   if     strncmpi(TS0,'st',2),  TS=0; %% 'standard'
   elseif strncmpi(TS0,'ha',2),  TS=1; %% 'harmonic'
   elseif strncmpi(TS0,'se',2),  TS=2; %% 'searchspace'
   elseif strncmpi(TS0,'bv',2),  TS=3;
   elseif strncmpi(TS0,'av',2),  TS=4;
   end
else
   TS=max(0,min(4,round(TS0(1,1))));
end

%------- set targets ----------------
if isempty(Sigma)
   if TS==1, Sigma=[0,1]; else, Sigma = 'LM'; end
elseif ischar(Sigma)
   switch Sigma
     case {'LM','LR','SR','BE','SM'}
     if ~ok, TS=3; end
   end
else
  [k,l]=size(Sigma);
  if l==1, Sigma=[Sigma,ones(k,1)]; l=2; end
  Sigma=ScaleEig(Sigma);
end

if ischar(Sigma) & TS<2 
   msg1=sprintf('   The choice sigma = ''%s'' does not match the\n',Sigma);
   msg2=sprintf('   selected test subspace. Specify a numerical\n');
   msg3=sprintf('   value for sigma (e.g. sigma = '); msg4='';
   switch Sigma
      case {'LM','LR'}
         msg4=sprintf('[1,0]');     
      case {'SM','SR','BE'}
         msg4=sprintf(' [0,1]');  
   end
   msg5=sprintf('),\n   or continue with ''TestSpace''=''B*V''.');
   msg=[msg1,msg2,msg3,msg4,msg5];
   button=questdlg(msg,'Targets and test subspaces','Continue','Stop','Continue');
   if strcmp(button,'Continue'), TS=3; else, n=-1; return, end
end

%------- linear solver --------------------------
lsolver  = findfield(options,fopts,'lso',lsolver);

method=strvcat('exact','olsen','iluexact','gmres','cgstab','bicgstab');
j=strmatch(lower(lsolver),method);
if isempty(j), 
  msg=['The linear solver ''',lsolver,''' is not included.'];
  msg=[msg,sprintf('\nIs GMRES ok?')];
  button=questdlg(msg,'Linear solver','Yes','No','Yes');
  if strcmp(button,'Yes'), j=4; ls_maxit=5; else, n=-1; return, end
end
if     j==1, lsolver='exact'; Precond_Form = 0;
elseif j==2 | j==3, lsolver='olsen';
elseif j==4, lsolver='gmres';  ls_maxit=5;
else,        lsolver='cgstab'; ls_tol=1.0e-10;
end

ls_maxit= findfield(options,fopts,'ls_m',ls_maxit);
ls_tol  = findfield(options,fopts,'ls_t',ls_tol);
ell     = findfield(options,fopts,'ls_e',ell);

par=[ls_tol,ls_maxit,ell];


%----- Display the parameters that are used ----------------
if DISP 

  fprintf('\n'),fprintf('PROBLEM\n')
  switch Operator_Form
     case {1,4,6,7} 
  fprintf('%13s: %s\n','A',StrOp(Operator_A));
  fprintf('%13s: %s\n','B',StrOp(Operator_B));
     case 2
  fprintf('%13s: ''%s''  ([Av,Bv] = %s(v))\n','A,B',Operator_A,Operator_A);
     case 3
  fprintf('%13s: ''%s''\n','A,B',Operator_A);
  fprintf('%15s(Av = %s(v,''A''), Bv = %s(v,''B''))\n',...
                                '',Operator_A,Operator_A);
     case {5,8}
  fprintf('%13s: %s\n','A',StrOp(Operator_A));
  fprintf('%13s: %s\n','B','Identity   (B*v = v)');
  end

  fprintf('%13s: %i\n','dimension',n);
  fprintf('%13s: %i\n\n','nselect',nselect);

  if length(jj)>0 & (ischar(Operator_A) | ischar(Operator_B))
    msgj=sprintf(', %i',jj); msgo=''; 
    if ischar(Operator_A), msgo=sprintf(' ''%s''',Operator_A); end
    if ischar(Operator_B), msgo=sprintf('%s ''%s''.',msgo,Operator_B); end
    fprintf(' The JDQZ input arguments, number%s, are\n',msgj)
    fprintf(' taken as input parameters 3:%i for%s.\n\n',length(jj)+2,msgo);
  end

  fprintf('TARGET\n')
  if ischar(Sigma)
     fprintf('%13s: ''%s''\n','sigma',Sigma)
  else 
     fprintf('%13s: %s\n','sigma',mydisp(Sigma(1,:)))
     for j=2:size(Sigma,1), 
        fprintf('%13s: %s\n','',mydisp(Sigma(j,:)))
     end
  end

  fprintf('\nOPTIONS\n')
  fprintf('%13s: %g\n','Schur',SCHUR)
  fprintf('%13s: %g\n','Tol',tol)
  fprintf('%13s: %i\n','Disp',DISP)
  fprintf('%13s: %i\n','jmin',jmin)
  fprintf('%13s: %i\n','jmax',jmax)
  fprintf('%13s: %i\n','MaxIt',maxit)
  fprintf('%13s: %s\n','v0',StrOp(V))
  fprintf('%13s: %i\n','Pairs',PAIRS)
  fprintf('%13s: %i\n','AvoidStag',JDV)
  fprintf('%13s: %i\n','NSigma',NSIGMA)
  fprintf('%13s: %g\n','Track',track)
  fprintf('%13s: %g\n','FixShift',FIX)
  fprintf('%13s: %i\n','Chord',CHORDALDISTANCE)
  fprintf('%13s: ''%s''\n','LSolver',lsolver)
  str=sprintf('%g ',ls_tol);
  fprintf('%13s: [ %s]\n','LS_Tol',str)
  fprintf('%13s: %i\n','LS_MaxIt',ls_maxit)
  if strcmp(lsolver,'cgstab')
    fprintf('%13s: %i\n','LS_ell',ell)
  end

  DisplayPreconditioner(n); fprintf('\n')
   
  switch TS
     case 0, str='Standard, W = alpha''*A*V + beta''*B*V';
     case 1, str='Harmonic, W = beta*A*V - alpha*B*V';
     case 2, str='SearchSpace, W = V';
     case 3, str='Petrov, W = B*V';
     case 4, str='Petrov, W = A*V';
   end
   fprintf('%13s: %s\n','TestSpace',str); fprintf('\n\n');

end

return
%------------------------------------------------------------------------
function msg=StrOp(Op)

  if ischar(Op)
    msg=sprintf('''%s''',Op);
  elseif issparse(Op), [n,k]=size(Op);
    msg=sprintf('[%ix%i sparse]',n,k);
  else, [n,k]=size(Op);
    msg=sprintf('[%ix%i double]',n,k);
  end

return
%------------------------------------------------------------------------
function DisplayPreconditioner(n)
  global Precond_Form Precond_Type ...
         Precond_L Precond_U Precond_P

  FP=Precond_Form;
  switch Precond_Form
    case 0,    
      fprintf('%13s: %s\n','Precond','No preconditioner'); 
    case {1,2,5}
      fprintf('%13s: %s','Precond',StrOp(Precond_L));
      if FP==2, fprintf('  (M\\v = %s(v,''preconditioner''))',Precond_L); end 
      fprintf('\n')
    case {3,4,6,7} 
      fprintf('%13s: %s\n','L precond',StrOp(Precond_L));
      fprintf('%13s: %s\n','U precond',StrOp(Precond_U));
      if FP==7, fprintf('%13s: %s\n','P precond',StrOp(Precond_P)); end
      if FP==4,fprintf('%15s(M\\v = %s(%s(v,''L''),''U''))\n',...
          '',Precond_L,Precond_L); end 
  end

  if FP
  switch Precond_Type
    case 0, str='explicit left';
    case 1, str='explicit right';
    case 2, str='implicit';
  end
  fprintf('%15sTo be used as %s preconditioner.\n','',str)
  end

return
%------------------------------------------------------------------------
function possibilities

 fprintf('\n')
 fprintf('PROBLEM\n')
 fprintf('            A: [ square matrix | string ]\n');
 fprintf('            B: [ square matrix {identity} | string ]\n');
 fprintf('      nselect: [ positive integer {5} ]\n\n');
 fprintf('TARGET\n')
 fprintf('        sigma: [ vector of scalars | \n');
 fprintf('                 pair of vectors of scalars |\n');
 fprintf('                 {''LM''} | {''SM''} | ''LR'' | ''SR'' | ''BE'' ]\n\n');

 fprintf('OPTIONS\n');
 fprintf('        Schur: [ yes | {no} ]\n');
 fprintf('          Tol: [ positive scalar {1e-8} ]\n');
 fprintf('         Disp: [ yes | {no} ]\n');
 fprintf('         jmin: [ positive integer {nselect+5} ]\n');
 fprintf('         jmax: [ positive integer {jmin+5} ]\n');
 fprintf('        MaxIt: [ positive integer {200} ]\n');
 fprintf('           v0: ');
 fprintf('[ size(A,1) by p vector of scalars {rand(size(A,1),1)} ]\n');
 fprintf('    TestSpace: ');
 fprintf('[ Standard | {Harmonic} | SearchSpace | BV | AV ]\n');
 fprintf('        Pairs: [ yes | {no} ] \n');
 fprintf('    AvoidStag: [ yes | {no} ]\n');
 fprintf('        Track: [ {yes} | no non-negative scalar {1e-4} ]\n');
 fprintf('       NSigma: [ yes | {no} ]\n');
 fprintf('     FixShift: [ yes | {no} | non-negative scalar {1e+3} ]\n');
 fprintf('        Chord: [ {yes} | no ]\n');
 fprintf('        Scale: [ positive scalar {1} ]\n');
 fprintf('      LSolver: [ {gmres} | bicgstab ]\n');
 fprintf('       LS_Tol: [ row of positive scalar {[0.7,0.49]} ]\n');
 fprintf('     LS_MaxIt: [ positive integer {5} ]\n');
 fprintf('       LS_ell: [ positive integer {4} ]\n');
 fprintf('      Precond: ');
 fprintf('[ n by n matrix {identity} | n by 2n matrix | string ]\n');
 fprintf('    L_Precond: same as ''Precond''\n');
 fprintf('    U_Precond: [ n by n matrix {identity} | string ]\n');
 fprintf('    P_Precond: [ n by n matrix {identity} ]\n');
 fprintf(' Type_Precond: [ {left} | right | implicit ]\n');    
 fprintf('\n')


return
%------------------------------------------------------------------------
function x = boolean(x,gamma,string)
%Y = BOOLEAN(X,GAMMA,STRING)
%  GAMMA(1) is the default. 
%  If GAMMA is not specified, GAMMA = 0.
%  STRING is a cell of accepted strings. 
%  If STRING is not specified STRING = {'n' 'y'}
%  STRING{I} and GAMMA(I) are accepted expressions for X 
%  If X=GAMMA(I) then Y=X. If the first L characters
%  of X matches those of STRING{I}, then Y=GAMMA(I+1).
%  Here L=SIZE({STRING{1},2).
%  For other values of X, Y=GAMMA(1);

if nargin < 2, gamma=[0,0,1]; end
if nargin < 3, string={'n' 'y'}; end

if ischar(x)
  l=size(string{1},2);
  i=min(find(strncmpi(x,string,l)));
  if isempty(i), i=0; end, x=gamma(i+1);
elseif max((gamma-x)==0)
elseif gamma(end) == inf
else, x=gamma(1);
end
  
return
%------------------------------------------------------------------------
function [a,ok]=findfield(options,fopts,str,default,gamma,stri)
% Searches the fieldnames in FOPTS for the string STR.
% The field is detected if only the first part of the fieldname
% matches the string STR. The search is case insensitive.
% If the field is detected, then OK=1 and A is the fieldvalue.
% Otherwise OK=0 and A=DEFAULT

   l=size(str,2); j=min(find(strncmpi(str,fopts,l)));
   if ~isempty(j)
      a=getfield(options,char(fopts(j,:))); ok=1;
      if nargin == 5, a = boolean(a,[default,gamma]); 
      elseif nargin == 6, a = boolean(a,[default,gamma],stri); end
   elseif nargin>3
      a=default; ok=0;
   else
      a=[]; ok=0;
   end

return
%%%======================================================================
function n=CheckMatrix(A,gamma)

  if ischar(A), n=-1;
    if exist(A) ~=2
      msg=sprintf('  Can not find the M-file ''%s.m''  ',A);
      errordlg(msg,'MATRIX'),n=-2;
    end
    if n==-1 & gamma, eval('n=feval(A,[],''dimension'');','n=-1;'), end
  else, [n,n] = size(A);
    if any(size(A) ~= n)
      msg=sprintf('  The operator must be a square matrix or a string.  ');
      errordlg(msg,'MATRIX'),n=-3;
    end
  end

return
%------------------------------------------------------------------------
function [Av,Bv,n]=CheckDimMV(v)
% [Av,Bv]=CheckDimMV(v)
%  Av=A*v, Bv=B*v Checks correctness operator definitions

global Operator_Form Operator_MVs Operator_A Operator_B Operator_Params

[n,k]=size(v);
if k>1, V=v(:,2:k); v=v(:,1); end

Bv=v;
switch Operator_Form
   case 1 % both Operator_A and B are strings 
      Operator_Params{1}=v;
      Av=feval(Operator_A,Operator_Params{:}); 
      Bv=feval(Operator_B,Operator_Params{:});
   case 2 %%% or Operator_Form=3 or Operator_Form=5???
      Operator_Params{1}=v;
      eval('[Av,Bv]=feval(Operator_A,Operator_Params{:});','Operator_Form=3;')
      if Operator_Form==3, 
        ok=1; Operator_Params{2}='A';
        eval('Av=feval(Operator_A,Operator_Params{:});','ok=0;')
        if ok, 
          Operator_Params{2}='B';
          eval('Bv=feval(Operator_A,Operator_Params{:});','ok=0;'), 
        end
        if ~ok, 
          Operator_Form=5; Operator_Params{2}=''; 
          Av=feval(Operator_A,Operator_Params{:}); 
        end 
      end
   case 3
      Operator_Params{1}=v;
      Operator_Params{2}='A';
      Av=feval(Operator_A,Operator_Params{:});
      Operator_Params{2}='B';
      Bv=feval(Operator_A,Operator_Params{:});
   case 4
      Operator_Params{1}=v;
      Av=feval(Operator_A,Operator_Params{:}); 
      Bv=Operator_B*v;
   case 5
      Operator_Params{1}=v;
      Av=feval(Operator_A,Operator_Params{:});
   case 6
      Av=Operator_A*v; 
      Operator_Params{1}=v;
      Bv=feval(Operator_B,Operator_Params{:});
   case 7
      Av=Operator_A*v; 
      Bv=Operator_B*v;
   case 8
      Av=Operator_A*v;
end

Operator_MVs = 1;

ok_A=min(size(Av)==size(v));
ok_B=min(size(Bv)==size(v));

if ok_A & ok_B
  if k>1,
    ok=1; eval('[AV,BV]=MV(V);','ok=0;') 
    if ok
      Av=[Av,AV]; Bv=[Bv,BV];
    else
      for j=2:k, [Av(:,j),Bv(:,j)]=MV(V(:,j-1)); end
    end 
  end
  return
end

if ~ok_A
  Operator=Operator_A;
elseif  ~ok_B
  Operator=Operator_B;
end  
msg=sprintf(' %s does not produce a vector of size %d',Operator,n)
errordlg(msg,'MATRIX'), n=-1; 

return
%------------------------------------------------------------------------
function n=SetPrecond(n,TP)
% finds out how the preconditioners are defined (Precond_Form)
% and checks consistency of the definitions.
%
% If M is the preconditioner then P*M=L*U. Defaults: L=U=P=I.
%
% Precond_Form
%       0:   no L
%       1:   L M-file, no U,     L ~= A
%       2:   L M-file, no U,     L == A
%       3:   L M-file, U M-file, L ~= A, U ~= A, L ~=U
%       4:   L M-file, U M-file, L == U
%       5:   L matrix, no U
%       6:   L matrix, U matrix  no P
%       7:   L matrix, U matrix, P matrix
%
% Precond_Type
%       0:   Explicit left
%       1:   Explicit right
%       2:   Implicit
%
  global Operator_A ...
         Precond_Form Precond_Solves ...
         Precond_Type ...
         Precond_L Precond_U Precond_P Precond_Params

  Precond_Type = 0;
  if ischar(TP)
    TP=lower(TP(1,1));
    switch TP
      case 'l'
        Precond_Type = 0;
      case 'r'
        Precond_Type = 1;
      case 'i'
        Precond_Type = 2; 
    end
  else
    Precond_Type=max(0,min(fix(TP),2));
  end

  Precond_Solves = 0;

  % Set type preconditioner
  Precond_Form=0;
  if isempty(Precond_L), return, end

  if ~isempty(Precond_U) & ischar(Precond_L)~=ischar(Precond_U)
    msg=sprintf('  L and U should both be strings or matrices');
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end
  if ~isempty(Precond_P) & (ischar(Precond_P) | ischar(Precond_L))
    msg=sprintf('  P can be specified only if P, L and U are matrices'); 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end  
  tp=1+4*~ischar(Precond_L)+2*~isempty(Precond_U)+~isempty(Precond_P);
  if tp==1, tp = tp + strcmp(Precond_L,Operator_A); end
  if tp==3, tp = tp + strcmp(Precond_L,Precond_U); end
  if tp==3 & strcmp(Precond_U,Operator_A)
    msg=sprintf('  If L and A use the same M-file,')
    msg=[msg,sprintf('\n  then so should U.')]; 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end
  if tp>5, tp=tp-1; end,  Precond_Form=tp;

  % Check consistency definitions
  if tp<5 & exist(Precond_L) ~=2
    msg=sprintf('  Can not find the M-file ''%s.m''  ',Precond_L); 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end

  ok=1; 
  if tp == 2
    Precond_Params{1}=zeros(n,1);
    Precond_Params{2}='preconditioner';
    eval('v=feval(Operator_A,Precond_Params{:});','ok=0;')
    if ~ok
       msg='Preconditioner and matrix use the same M-file';
       msg=[msg,sprintf(' ''%s.m''  \n',Precond_L)];
       msg=[msg,'Therefore the preconditioner is called as'];
       msg=[msg,sprintf('\n\n\tw=%s(v,''preconditioner'')\n\n',Precond_L)]; 
       msg=[msg,sprintf('Put this "switch" in ''%s.m''.',Precond_L)];
    end
  end

  if tp == 4 | ~ok
    ok1=1;
    Precond_Params{1}=zeros(n,1);
    Precond_Params{2}='L';
    eval('v=feval(Precond_L,Precond_Params{:});','ok1=0;')
    Precond_Params{2}='U';
    eval('v=feval(Precond_L,Precond_Params{:});','ok1=0;')
    if ok1 
      Precond_Form = 4; Precond_U = Precond_L; ok=1;
    else
      if tp == 4
        msg='L and U use the same M-file';
        msg=[msg,sprintf(' ''%s.m''   \n',Precond_L)];
        msg=[msg,'Therefore L and U are called'];
        msg=[msg,sprintf(' as\n\n\tw=%s(v,''L'')',Precond_L)]; 
        msg=[msg,sprintf(' \n\tw=%s(v,''U'')\n\n',Precond_L)]; 
        msg=[msg,sprintf('Check the dimensions and/or\n')];
        msg=[msg,sprintf('put this "switch" in ''%s.m''.',Precond_L)];
      end
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    end
  end

  if tp==1 | tp==3
    Precond_Params{1}=zeros(n,1); Precond_Params{2}='';
    eval('v=feval(Precond_L,Precond_Params{:});','ok=0')
    if ~ok
       msg=sprintf('''%s'' should produce %i-vectors',Precond_L,n); 
       errordlg(msg,'PRECONDITIONER'), n=-1; return 
    end
  end

  if tp==3
    if exist(Precond_U) ~=2
      msg=sprintf('  Can not find the M-file ''%s.m''  ',Precond_U);
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    else
      Precond_Params{1}=zeros(n,1); Precond_Params{2}='';
      eval('v=feval(Precond_U,Precond_Params{:});','ok=0')
      if ~ok
        msg=sprintf('''%s'' should produce %i-vectors',Precond_U,n);
        errordlg(msg,'PRECONDITIONER'), n=-1; return 
      end
    end
  end
  
  if tp==5
    if min([n,2*n]==size(Precond_L)) 
      Precond_U=Precond_L(:,n+1:2*n); Precond_L=Precond_L(:,1:n); 
      Precond_Form=6;
    elseif min([n,3*n]==size(Precond_L)) 
      Precond_U=Precond_L(:,n+1:2*n); Precond_P=Precond_L(:,2*n+1:3*n);
      Precond_L=Precond_L(:,1:n); Precond_Form=7;
    elseif ~min([n,n]==size(Precond_L)) 
      msg=sprintf('The preconditioning matrix\n');
      msg2=sprintf('should be %iX%i or %ix%i ([L,U])\n',n,n,n,2*n); 
      msg3=sprintf('or %ix%i ([L,U,P])\n',n,3*n); 
      errordlg([msg,msg2,msg3],'PRECONDITIONER'), n=-1; return
    end
  end

  if tp==6 & ~min([n,n]==size(Precond_L) & [n,n]==size(Precond_U))
    msg=sprintf('Both L and U should be %iX%i.',n,n); n=-1;
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

  if tp==7 & ~min([n,n]==size(Precond_L) & ...
       [n,n]==size(Precond_U) & [n,n]==size(Precond_P))
    msg=sprintf('L, U, and P should all be %iX%i.',n,n); n=-1;
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

return
%%%======================================================================
%%%========= DISPLAY FUNCTIONS ===========================================
%%%======================================================================
function Result(Sigm,Target,S,T,tol)

  if nargin == 1
     fprintf('\n\n %20s\n','Eigenvalues')
     for j=1:size(Sigm,1), 
       fprintf('%2i %s\n',j,mydisp(Sigm(j,:),5))
     end
     return
  end

  fprintf('\n\n%17s  %25s%18s\n','Targets','Eigenvalues','RelErr/Cond')
  for j=1:size(S,1), 
    l=Target(j,1); s=S(j,:); t=T(j,:); lm=mydisp(s/t,5);
    if l>0
      fprintf('%2i %8s ''%s'' %9s %23s',j,'',Sigm(l,:),'',lm)
    else
      fprintf('%2i %23s %23s',j,mydisp(Target(j,2:3),5),lm)
    end
    re=1-tol/abs(t); 
    if re>0, re=(1+tol/abs(s))/re; re=min(re-1,1); else, re=1; end
    if re>0.999, fprintf('  1\n'), else, fprintf('  %5.0e\n',re), end
  end

return

%------------------------------------------------------------------------
function s=mydisp(lambda,d)

  if nargin <2, d=4; end

  if size(lambda,2)==2,
    if lambda(:,2)==0, s='Inf'; return, end
    lambda=lambda(:,1)/lambda(:,2); 
  end

  a=real(lambda); b=imag(lambda);
  if max(a,b)==inf, s='Inf'; return, end

  e=0; 
  if abs(a)>0; e= floor(log10(abs(a))); end
  if abs(b)>0; e= max(e,floor(log10(abs(b)))); end
  p=10^d; q=p/(10^e); a=round(a*q)/p; b=round(b*q)/p;

  st=['%',num2str(d+2),'.',num2str(d),'f'];
  ab=abs(b)==0;
  if ab, str=[''' ']; else, str=['''(']; end
  if a>=0, str=[str,'+']; a=abs(a); end, str=[str,st];
  if ab, str=[str,'e']; else
  if b>=0, str=[str,'+']; b=abs(b); end, 
  str=[str,st,'i)e']; end
  if e>=0, str=[str,'+']; else, str=[str,'-']; e=-e; end
  if e<10, str=[str,'0%i''']; else, str=[str,'%i''']; end
  if ab, s=eval(sprintf(str,a,e));
  else, s=eval(sprintf(str,a,b,e)); end

return
%%%======================================================================
function ShowEig(theta,target,k)

  if ischar(target)
    fprintf('\n target=''%s'', ',target)
  else
    fprintf('\n target=%s, ',mydisp(target,5))
  end
  fprintf('lambda_%i=%s\n',k,mydisp(theta,5)); 

return
%%%======================================================================
function texttest(s,nr,gamma)

  if nargin<3, gamma=100*eps; end

  if nr > gamma
 % if nr>100*eps
     fprintf('\n %35s is: %9.3g\t',s,nr)
  end

return
%%%======================================================================


