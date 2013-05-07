function [ Z, W, E ] = l1_low_rank( X, beta, lambda, maxIter )
%L1_LOW_RANK Summary of this function goes here
%   Detailed explanation goes here

[d,n]=size(X);
m=n;
fprintf(1,'n is %d\n',n);
tol1 = 1e-4;%threshold for the error in constraint
tol2 = 1e-5;%threshold for the change in the solutions


Z=zeros(m,n);
W=Z;
E=zeros(d,n);
Y1=zeros(d,n);
Y2=zeros(m,n);

mu = min(d,n)*tol2;
max_mu=1e10;
rho_0=1.1;
rho=rho_0;

epsilon1=1e-6;
epsilon2=1e-2;

norm2X = norm(X,2);
A=X;
eta1 = norm2X*norm2X*1.02;%eta needs to be larger than ||X||_2^2, but need not be too large.
fprintf(1,'eta1 is %f\n',eta1);
k=0;

Zk_1=Z;
Wk_1=W;
Ek_1=E;

Xf=norm(X,'fro');
AZ=zeros(d,n);

MAX_ITER=maxIter;
iter=0;

sv = 5;

% main loop
while true
    tic;
    % update Z
    Zk_1=Z;
    T=Z+(A'*(X-AZ-E+Y1/mu)-(Z-W+Y2/mu))/eta1;
    % if choosvd(n, sv) == 1
        % % fprintf(1,'partial svd use lansvd\n');
        % [U S V] = lansvd(T, sv, 'L');
    % else
        % % fprintf(1,'matlab svd\n');
        % [U S V] = svd(T, 'econ');
    % end
    % diagS = diag(S);
    % svp = length(find(diagS > 1/(eta1*mu)));
    % if svp < sv
        % sv = min(svp + 1, n);
    % else
        % sv = min(svp + round(0.05*n), n);
    % end
    % Z = U(:, 1:svp) * diag(diagS(1:svp) - 1/(mu*eta1)) * V(:, 1:svp)';

    [Z,svp]=singular_value_shrinkage(T,1/(eta1*mu));
    Z=max(Z,0); % convex 
    Z=Z-diag(diag(Z)); % TODO: not tested, added 2013-03-28 diag(Z)=0

    AZ=A*Z;
    % update W
    Wk_1=W;
    W=max(W,0);
    W=wthresh(Z+Y2/mu,'s',beta*(1/mu)); 
    W=W-diag(diag(W)); % TODO: not tested, added 2013-03-28 diag(W)=0
    % W=max(wthresh(Z+Y2/mu,'s',beta*mu^-1),0); 
    % update E
    Ek_1=E;
    E=l21(X-AZ+Y1/mu,lambda*(1/mu));
    
    % update lagrange multipliers
    Y1=Y1+mu*(X-AZ-E);
    Y2=Y2+mu*(Z-W);
    
    % stop criteria
    t1=norm(X-AZ-E,'fro')/Xf;
    t2=mu*max(max(sqrt(eta1)*norm(Z-Zk_1,'fro'),norm(W-Wk_1,'fro')),norm(E-Ek_1,'fro'))/Xf;

    % update mu
    if t2>=epsilon2
        rho=1;
    else
        rho=rho_0;
    end
    mu = min(max_mu,rho*mu);
    
    %update k
    k=k+1;
    t=toc;
    % fprintf(1,'one iteration takes:%f\n',t);
    
    iter=iter+1;
    if iter==1 || mod(iter,50)==0 
        fprintf(1,'iter %d,svp %d,reconstruction error is %f,change error is %f\n',iter,svp,t1,t2);
    end

    % 判断中止条件
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        break;
    end
    
    % terminate condition
    if t1>=epsilon1 || t2>=epsilon2
    else
        fprintf(1,'stop criteria reached! terminate!\n');
        break;
    end
end

