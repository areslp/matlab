function [ Z, E ] = low_rank( X, lambda, maxIter )
%LOW_RANK Summary of this function goes here
%   Detailed explanation goes here

d=size(X,1);
n=size(X,2);

Z=zeros(n,n);
J=Z;

E=zeros(d,n);
Y1=zeros(d,n);
Y2=zeros(n,n);

mu=1e-6;
max_mu=10^10;
rho=1.1;
epsilon=1e-8;

xtx=X'*X;
inv_x=inv(xtx+eye(n));

MAX_ITER=maxIter;
iter=0;

while true
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        break;
    end
    % 1. update J
    tic;
    Y=Z+Y2/mu;
    tau=1/mu;
    J=singular_value_shrinkage(Y,tau);
    t=toc;
    fprintf(1,'step 1: update J takes: %f\n',t);
    % 2. update Z
    tic;
    Z=inv_x*(xtx-X'*E+J+(X'*Y1-Y2)/mu);
    t=toc;
    fprintf(1,'step 2: update Z takes: %f\n',t);
    % 3. update E
    tic;
    tlambda=lambda/mu;
    Q=X-X*Z+Y1/mu;
    E=l21(Q,tlambda);
    t=toc;
    fprintf(1,'step 3: update E takes: %f\n',t);
    % 4. update Y1, Y2
    tic;
    xz=X-X*Z;
    zj=Z-J;
    leq1=xz-E;
    leq2=zj;
    Y1=Y1+mu*(leq1);
    Y2=Y2+mu*(leq2);
    t=toc;
    fprintf(1,'step 4: update Y1, Y2 takes: %f\n',t);
    % 5. update mu
    tic;
    mu=min(rho*mu,max_mu);
    t=toc;
    fprintf(1,'step 5: update mu takes: %f\n',t);
    % 6. check the convergence
    if max(max(abs(leq1)))<epsilon && max(max(abs(leq2)))<epsilon
        break;
    end
    iter=iter+1;
end

end

