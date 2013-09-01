function [Z] = SSQP(X, lambda)
% implement the algorithm in paper "Efficient Subspace Segmentation via Quadratic Programming"
% f(Z)=||XZ-X||_F^2+lambda e^TZ^TZe

[m,n]=size(X);
Z=rand(n,n);
Zk=Z;
e=ones(n,1);
E=ones(n,n);
tau=0.0;
rho=0.0;
XtX=X'*X;
tol1=1e-4;
tol2=1e-1;

convergenced=false;
iter=0;
while ~convergenced
    gz=gradient(XtX,Z,E,lambda);
    D=projectZ(Z-tau*gz)-Z;
    % compute rho by line search
    % tic
    % rho=backTracking(Z,D,X,XtX,lambda,E,@gradient,@fx);
    rho=0.5;
    % t=toc;
    % fprintf(1,'back tracking takes: %f, rho is %f\n',t,rho);
    Zk=Z;
    Z=Z+rho*D;
    s=reshape(Z-Zk,n*n,1);
    gzn=gradient(XtX,Z,E,lambda);
    y=reshape(gzn-gz,n*n,1);
    % y'*y
    % s'*y
    tau=(s'*s)/(s'*y);
    % check convergence
    cc1=norm(Z-Zk,'fro');
    cc2=tau;
    if cc1<tol1 && cc2<tol2
        convergenced=true;
    end
    if mod(iter,100)==0 || convergenced
        fprintf(1,'iter is %d, cc1 is %f, cc2 is %f\n',iter,cc1,cc2);
    end
    iter=iter+1;
end


function [Z] = projectZ(Z)
    Z=max(Z,0);

function [g] = gradient(XtX,Z,E,lambda)
    g=2*XtX*Z-2*XtX+2*lambda*Z*E;

function [f] = fx(X,Z,lambda)
    f=norm(X*Z-X,'fro')+lambda*norm(Z'*Z,1);

function [alpha] = backTracking(Z,dir,X,XtX,lambda,E,gradient,fx)
    c1=1e-4;
    c2=0.1;
    tau=0.5;
    alpha_min=1e-8;
    % init guess of alpha
    alpha=1.0;
    k=0;
    convergenced=false;
    while ~convergenced
        % check Armijo rule and curvature
        Zn=Z+alpha*dir; 
        % fx
        fz=fx(X,Z,lambda);
        fzn=fx(X,Zn,lambda);
        % gx
        gz=gradient(XtX,Z,E,lambda);
        gzn=gradient(XtX,Zn,E,lambda);
        cond1=fzn-fz-c1*alpha*dir'*gz;
        cond2=dir'*gzn-c2*dir'*gz;
        if (max(cond1(:))<=0 && min(cond2(:))>=0)
            convergenced=true;
        end
        if alpha<alpha_min
            convergenced=true;
        end
        % update
        alpha=alpha*tau;
        k=k+1;
    end

