function x = lsqrSOLtest( m, n, damp )

%      x = lsqrSOLtest(  m, n, damp);
%      x = lsqrSOLtest( 10,10, 0   );
%      x = lsqrSOLtest( 20,10, 0   );
%      x = lsqrSOLtest( 20,10, 0.1 );
%
% If m = n and damp = 0, this sets up a system Ax = b
% and calls lsqrSOL.m to solve it.  Otherwise, the usual
% least-squares or damped least-squares problem is solved.

% 11 Apr 1996: First version for distribution with lsqr.m.
% 07 Aug 2002: LSQR's output parameter rnorm changed to r1norm, r2norm.
% 03 May 2007: Allow A to be a matrix or a function handle.
%              Private function Aprodxxx defines matrix-vector products
%              for a specific A.
% 24 Dec 2010: A*v and A'*v use inputs (v,1) and (v,2), not (1,v) and (2,v).

%              Michael Saunders, Systems Optimization Laboratory,
%              Dept of MS&E, Stanford University.
%-----------------------------------------------------------------------

  A      = @(v,mode) Aprodxxx( v,mode,m,n );  % Nested function

  xtrue  = (n : -1: 1)';
  b      = A(xtrue,1);

  atol   = 1.0e-6;
  btol   = 1.0e-6;
  conlim = 1.0e+10;
  itnlim = 10*n;
  show   = 1;

  [ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ] ...
      = lsqrSOL( m, n, A, b, damp, atol, btol, conlim, itnlim, show );

  disp(' ');   j1 = min(n,5);   j2 = max(n-4,1);
  disp('First elements of x:');  disp(x(1:j1)');
  disp('Last  elements of x:');  disp(x(j2:n)');

  r    = b - A(x,1);
  r1   = norm(r);
  r2   = norm([r; (-damp*x)]);
  disp(' ')
  str1 = sprintf( 'r1norm, r2norm (est.)  %10.3e %10.3e', r1norm, r2norm );
  str2 = sprintf( 'r1norm, r2norm (true)  %10.3e %10.3e', r1    , r2     );
  disp(str1)
  disp(str2)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Nested functions (only 1 here).
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function y = Aprodxxx( x, mode, m, n )

    % Private function.
    % if mode = 1, computes y = A*x
    % if mode = 2, computes y = A'*x
    % for some matrix  A.
    %
    % This is a simple example for testing  LSQR.
    % It uses the leading m*n submatrix from
    % A = [ 1
    %       1 2
    %         2 3
    %           3 4
    %             ...
    %               n ]
    % suitably padded by zeros.
    %
    % 11 Apr 1996: First version for distribution with lsqr.m.
    %              Michael Saunders, Dept of EESOR, Stanford University.

    if mode == 1
      d  = (1:n)';  % Column vector
      y1 = [d.*x; 0] + [0;d.*x];
      if m <= n+1
	y = y1(1:m);
      else         
	y = [     y1; 
	  zeros(m-n-1,1)];
      end
    else
      d  = (1:m)';  % Column vector
      y1 = [d.*x] + [d(1:m-1).*x(2:m); 0];
      if m >= n
	y = y1(1:n);
      else
	y = [y1;
	  zeros(n-m,1)];
      end
    end

  end % nested function Aprodxxx

end % function lsqrSOLtest


