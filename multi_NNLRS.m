function [Z,ZZ,E] = multi_NNLRS(X,lambda,beta,alpha)
% solve \sum_{i=1}^k(||Z_i||_*+beta||Z_i||_1+lambda||E_i||_{2,1})+alpha||ZZ||_{2,1}

tic;

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
    fprintf(1,'eta1{%d} is %f\n',i,eta1{i});
end
mu=1e-6;
max_mu=10^10;
rho=1.9;
% epsilon=1e-4;
% epsilon2=1e-5; % must be small!
epsilon=1e-6;
epsilon2=1e-2; % must be small!
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
	% for i=1:k
		% fprintf(1,'S{%d}, max: %f, min: %f\n',i,max(max(S{i})),min(min(S{i})));
	% end
    % update J_i
    Jk=J;
    [J]=cellfun(@updateJ,Z,J,Y2,cmu,cbeta,'UniformOutput',false);
    % for i=1:k
        % norm(Jk{i}-J{i})
    % end
    % update Z
    [F]=cellfun(@updateF,J,Y2,S,Y3,cmu,'UniformOutput',false); 

    % normalize matrix before L21, then restore them
    % CO=F;
    % for i=1:k
        % FN=sqrt(sum(F{i}.^2,1));
        % CO{i}=FN; % CO is the column norm of matrix F
        % F{i}=mnormalize_col(F{i});
    % end
    % save_matrix;
        
    [M]=cellfun(@updateM,F,'UniformOutput',false);
    MM=zeros(k,n*n);
    for i=1:k
		% TODO: normalize
		% fprintf(1,'M{%d}, max: %f, min: %f\n',i,max(max(M{i})),min(min(M{i})));
        % M{i} = (M{i} - min(M{i}(:))) ./ (max(M{i}(:))-min(M{i}(:)));
		% fprintf(1,'M{%d},max: %f, min: %f\n',i,max(max(M{i})),min(min(M{i})));
        MM(i,:)=M{i};
    end
    ZZ=l21(MM,alpha/(2*mu));
    % if alpha==0
        % assert(nnz(ZZ-MM)==0);
    % end
    % update Z_i
    Zk=Z;
    for i=1:k
        Z{i}=reshape(ZZ(i,:),n,n)';
        % Z{i}=Z{i}-diag(diag(Z{i}));
        % Z{i}=max(Z{i},0);

        % multiple the CO matrix to Z's columns
        % Z{i}=Z{i}*repmat(CO{i},size(Z{i},1),1);
    end
    % update E_i
    Ek=E;
    [E]=cellfun(@updateE,X,S,E,Y1,cmu,clambda,'UniformOutput',false);

    % parameter update rule

    % check convergence
    [Xv,Xc,ZJv,ZJc,ZSv,ZSc,Zc,Jc,Sc,Ec,Cmax] = cellfun(@caculateTempVars,X,S,E,Z,J,Zk,Jk,Sk,Ek,Xf,eta1,cmu,'UniformOutput',false);
    changeX=max([Xv{:}]);
    changeZJ=max([ZJv{:}]);
    changeZS=max([ZSv{:}]);
    gap=max([Cmax{:}]);
    if mod(iter,50)==0
        fprintf(1,'===========================================================================================================\n');
        fprintf(1,'gap between two iteration is %e,mu is %e\n',gap,mu);
        fprintf(1,'iter %d,mu is %e,ResidualX is %e,changeZJ is %e,changeZS is %e\n',iter,mu,changeX,changeZJ,changeZS);
        for i=1:k
            fprintf(1,'svp%d %d,',i,svp{i});
        end
        fprintf(1,'\n');
    end
    if changeX <= epsilon && gap <=epsilon2
    % if changeX <= epsilon && gap <=epsilon2 && changeZJ <= epsilon && changeZS <= epsilon
        convergenced=true;
        fprintf(2,'convergenced, iter is %d\n',iter);
        fprintf(2,'iter %d,mu is %e,ResidualX is %e,changeZJ is %e,changeZS is %e\n',iter,mu,changeX,changeZJ,changeZS);
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
    % save_matrix(J,S,Z,iter);
    iter=iter+1;
end

toc;

function [S,svp] = updateS(xtx,X,E,Y1,Z,S,Y3,eta1,mu)
    T=-mu*(xtx-xtx*S-X'*E+X'*Y1/mu+Z-S+Y3/mu);
    % argmin_{S} 1/(mu*eta1)||S||_*+1/2*||S-S_k+T/(mu*eta1)||_F^2
    [S,svp]=singular_value_shrinkage(S-T/(mu*eta1),1/(mu*eta1)); % TODO: sometimes PROPACK is slower than full svd, and sometimes it will throw the following error
    % S=S-diag(diag(S));
    % S=max(S,0);

function [J] = updateJ(Z,J,Y2,mu,beta)
    J=wthresh(Z+Y2/mu,'s',beta/mu); 
    % J=J-diag(diag(J));
    % J=max(J,0);

function [F] = updateF(J,Y2,S,Y3,mu)
    F=0.5*(J-Y2/mu+S-Y3/mu);

function [M] = updateM(F)
    n=length(F);
    M=reshape(F',1,n*n);

function [E] = updateE(X,S,E,Y1,mu,lambda)
    E=l21(X-X*S+Y1/mu,lambda/mu); % TODO: -E not E

function [Xv,Xc,ZJv,ZJc,ZSv,ZSc,Zc,Jc,Sc,Ec,Cmax] = caculateTempVars(X,S,E,Z,J,Zk,Jk,Sk,Ek,Xf,eta1,mu)
    Xc=X-X*S-E; 
    ZJc=Z-J;
    ZSc=Z-S;
    Xv=norm(Xc,'fro')/Xf;
    ZJv=norm(ZJc,'fro')/Xf;
    ZSv=norm(ZSc,'fro')/Xf;

    Zc=norm(Z-Zk,'fro')/Xf;
    Jc=norm(J-Jk,'fro')/Xf;
    Sc=norm(S-Sk,'fro')/Xf;
    Ec=norm(E-Ek,'fro')/Xf;

    Cmax=mu*(max([sqrt(eta1)*Sc Jc Zc Ec]));


function [Y1] = updateY1(Y1,mu,Xc)
    Y1=Y1+mu*Xc;
    
function [Y2] = updateY2(Y2,mu,ZJc)
    Y2=Y2+mu*ZJc;

function [Y3] = updateY3(Y3,mu,ZSc)
    Y3=Y3+mu*ZSc;


function [] = save_matrix(J,S,Z,iter)
    save(['m' num2str(iter) '.mat'],'J','S','Z');
    k=length(J);
    for i=1:k

    end
