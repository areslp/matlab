function [Z,ZZ,E] = multi_NNLRS(X,lambda,beta,alpha)
% init vars
k=length(X);
[m,n]=size(X{1});

Z=cell(k,1);
Z(1:k)={zeros(n)};
E=cell(k,1);
E(1:k)={zeros(m,n)};
S=cell(k,1);
S(1:k)={zeros(n)};
J=cell(k,1);
J(1:k)={zeros(n)};
Y1=cell(k,1);
Y1(1:k)={zeros(m,n)};
Y2=cell(k,1);
Y2(1:k)={zeros(n)};
Y3=cell(k,1);
Y3(1:k)={zeros(n)};
Zk=Z;
Ek=E;
Sk=S;
Jk=J;
svp=cell(k,1);
svp(1:k)={0};
F=Z;
ZZ=zeros(k,n*n);

% precomputed values
xtx=cell(k,1);
for i=1:k
    xtx{i}=X{i}'*X{i};
end
invx=cell(k,1);
for i=1:k
    invx{i}=inv(xtx{i}+eye(n));
end
Xf=cell(k,1);
for i=1:k
    Xf{i}=norm(X{i},'fro');
end
% the residual error and the error between Z,J,S
Xc=cell(k,1);
ZJc=cell(k,1);
ZSc=cell(k,1);

% parameters
norm2X=cell(k,1);
for i=1:k
    norm2X{i}=norm(X{i},2);
end
eta1=cell(k,1);
for i=1:k
    eta1{i}=norm2X{i}*norm2X{i}*1.02;%eta needs to be larger than ||X||_2^2, but need not be too large.
end
mu=1e-6;
max_mu=10^10;
rho=1.9;
% epsilon=1e-4;
% epsilon2=1e-5; % must be small!
epsilon=1e-6;
epsilon2=1e-5; % must be small!
MAX_ITER=1000;
iter=0;
convergenced=false;
clambda=cell(k,1);
clambda(1:k)={lambda};
cbeta=cell(k,1);
cbeta(1:k)={beta};

while ~convergenced
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        break;
    end
    cmu=cell(k,1);
    cmu(1:k)={mu};
    % update S_i
    Sk=S;
    [S, svp]=cellfun(@updateS,xtx,X,E,Y1,Z,S,Y3,eta1,cmu,'UniformOutput',false);
    % update J_i
    Jk=J;
    [J]=cellfun(@updateJ,Z,J,Y2,cmu,cbeta,'UniformOutput',false);
    % update Z
    [F]=cellfun(@updateF,J,Y2,S,Y3,cmu,'UniformOutput',false);
    [M]=cellfun(@updateM,F,'UniformOutput',false);
    for i=1:k
        ZZ(i,:)=M{i};
    end
    ZZ=l21(ZZ,alpha/(2*mu));
    % update Z_i
    Zk=Z;
    for i=1:k
        Z{i}=reshape(ZZ(i,:),n,n)';
    end
    % update E_i
    [E]=cellfun(@updateE,X,S,E,Y1,cmu,clambda,'UniformOutput',false);

    % parameter update rule

    % check convergence
    [Xv,Xc,ZJv,ZJc,ZSv,ZSc,Zc,Jc,Sc,Ec] = cellfun(@caculateTempVars,X,S,E,Z,J,Zk,Jk,Sk,Ek,Xf,'UniformOutput',false);
    changeX=max([Xv{:}]);
    changeZJ=max([ZJv{:}]);
    changeZS=max([ZSv{:}]);
    changeZ=max([Zc{:}]); 
    changeJ=max([Jc{:}]);
    changeS=max([Sc{:}]);
    changeE=max([Ec{:}]);
    tmp=[changeZ changeJ changeS changeE ];
    gap=mu*max(tmp);
    if mod(iter,50)==0
        fprintf(1,'===========================================================================================================\n');
        fprintf(1,'gap between two iteration is %f,mu is %f\n',gap,mu);
        fprintf(1,'iter %d,mu is %f,ResidualX is %f,changeZJ is %f,changeZS is %f\n',iter,mu,changeX,changeZJ,changeZS);
        for i=1:k
            fprintf(1,'svp%d %d,',i,svp{i});
        end
        fprintf(1,'\n');
    end
    % if changeX <= epsilon && changeZJ <= epsilon && changeZS <= epsilon
    if changeX <= epsilon && gap <=epsilon2 && changeZJ <= epsilon && changeZS <= epsilon
        convergenced=true;
        fprintf(2,'convergenced, iter is %d\n',iter);
        fprintf(2,'iter %d,mu is %f,ResidualX is %f,changeZJ is %f,changeZS is %f\n',iter,mu,changeX,changeZJ,changeZS);
        for i=1:k
            fprintf(1,'svp%d %d,',i,svp{i});
        end
        fprintf(1,'\n');
    end
    % update multipliers
    [Y1]=cellfun(@updateY1,Y1,cmu,Xc,'UniformOutput',false);
    [Y2]=cellfun(@updateY2,Y2,cmu,ZJc,'UniformOutput',false);
    [Y3]=cellfun(@updateY3,Y3,cmu,ZSc,'UniformOutput',false);
    % update parameters
    if gap < epsilon2
        mu=min(rho*mu,max_mu);
    end
    iter=iter+1;
end

function [S,svp] = updateS(xtx,X,E,Y1,Z,S,Y3,eta1,mu)
    T=-mu*(xtx-xtx*S-X'*E+X'*Y1/mu+Z-S+Y3/mu);
    % argmin_{S} 1/(mu*eta1)||S||_*+1/2*||S-S_k+T/(mu*eta1)||_F^2
    [S,svp]=singular_value_shrinkage(S-T/(mu*eta1),1/(mu*eta1)); % TODO: sometimes PROPACK is slower than full svd, and sometimes it will throw the following error

function [J] = updateJ(Z,J,Y2,mu,beta)
    J=wthresh(Z+Y2/mu,'s',beta/mu); 

function [RET] = updateF(J,Y2,S,Y3,mu)
    RET=1/2*(J-Y2/mu+S-Y3/mu);

function [M] = updateM(F)
    n=length(F);
    M=reshape(F',1,n*n);

function [E] = updateE(X,S,E,Y1,mu,lambda)
    E=l21(X*S-X-Y1/mu,lambda/mu);

function [Xv,Xc,ZJv,ZJc,ZSv,ZSc,Zc,Jc,Sc,Ec] = caculateTempVars(X,S,E,Z,J,Zk,Jk,Sk,Ek,Xf)
    Xc=X-X*S-E; 
    ZJc=Z-J;
    ZSc=Z-S;
    Xv=norm(Xc,'fro')/Xf;
    ZJv=norm(ZJc,'fro')/Xf;
    ZSv=norm(ZSc,'fro')/Xf;

    Zc=norm(Zk-Z,'fro')/Xf;
    Jc=norm(Jk-J,'fro')/Xf;
    Sc=norm(Sk-S,'fro')/Xf;
    Ec=norm(Ek-E,'fro')/Xf;

function [Y1] = updateY1(Y1,mu,Xc)
    Y1=Y1+mu*Xc;
    
function [Y2] = updateY2(Y2,mu,ZJc)
    Y2=Y2+mu*ZJc;

function [Y3] = updateY3(Y3,mu,ZSc)
    Y3=Y3+mu*ZSc;
