function [ZZ,Z,E] = multi_lrr(X,lambda,alpha)
% implement the algorithm described in paper "Multi-task Low-rank Affinity Pursuit for Image Segmentation"
% solve \sum_{i=1}^k(||J_i||_*+lambda||E_i||_{2,1})+alpha||ZZ||_{2,1} s.t. X=XS+E,Z=J,Z=S

k=size(X,1); % k*1 cell array, k views
[m,n]=size(X{1}); % every view has the same dimension
% initial matrix cell array
E=cell(k,1);
for i=1:k
    E{i}=zeros(m,n);
end
J=cell(k,1);
for i=1:k
    J{i}=zeros(n,n);
end
S=cell(k,1);
for i=1:k
    S{i}=zeros(n,n);
end
Z=cell(k,1);
for i=1:k
    Z{i}=zeros(n,n);
end
W=cell(k,1);
for i=1:k
    W{i}=zeros(n,n);
end
Y=cell(k,1);
for i=1:k
    Y{i}=zeros(m,n);
end
V=cell(k,1);
for i=1:k
    V{i}=zeros(n,n);
end
ZZ=zeros(k,n*n);
% k's iteration vars
Ek=cell(k,1);
for i=1:k
    Ek{i}=zeros(m,n);
end

Jk=cell(k,1);
for i=1:k
    Jk{i}=zeros(n,n);
end

Sk=cell(k,1);
for i=1:k
    Sk{i}=zeros(n,n);
end

Zk=cell(k,1);
for i=1:k
    Zk{i}=zeros(n,n);
end

% parameters
mu=1e-6;
max_mu=10^10;
rho=1.9;
% epsilon=1e-4;
% epsilon2=1e-5; % must be small!
epsilon=1e-8; % 1e-5
epsilon2=1e-8; % must be small! 1e-4

% pre caculate matrix value
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


sv=cell(k,1);
for i=1:k
    sv{i}=0;
end

svp=cell(k,1);
for i=1:k
    svp{i}=0;
end
F=cell(k,1);


M=cell(k,1);

MAX_ITER=100000000;
iter=0;
convergenced=false;


clambda=cell(k,1);
clambda(1:k)={lambda};

tic
while ~convergenced
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        save_matrix(J,S,Z,F);
        break;
    end
    cmu=cell(k,1);
    cmu(1:k)={mu};
    % update J_i
    Jk=J;
    [J, svp, sv]=cellfun(@updateJ,Z,W,cmu,sv,'UniformOutput',false);
    % update S_i
    Sk=S;
    S=cellfun(@updateS,invx,xtx,X,E,Z,Y,V,W,cmu,'UniformOutput',false,'ErrorHandler',@errorfun);
    % update ZZ
    [F]=cellfun(@updateF,J,S,W,V,cmu,'UniformOutput',false);
    [M]=cellfun(@updateM,F,'UniformOutput',false);
    MM=zeros(k,n*n);
    for i=1:k
		% TODO: normalize
		% fprintf(1,'M{%d}, max: %f, min: %f\n',i,max(max(M{i})),min(min(M{i})));
        % M{i} = (M{i} - min(M{i}(:))) ./ (max(M{i}(:))-min(M{i}(:)));
		% fprintf(1,'M{%d},max: %f, min: %f\n',i,max(max(M{i})),min(min(M{i})));
		MM(i,:)=M{i};
    end
	% fprintf(1,'============================================================\n');
    ZZ=l21(MM,alpha/(2*mu));
    % update Z_i
    for i=1:k
        Zk{i}=Z{i};
        Z{i}=reshape(ZZ(i,:),n,n)';
        % fprintf(1,'Z changed by L21, the diff between Z%d{%d} and Z{%d} is %e\n',iter,i,i,norm(Zk{i}-Z{i}));
        % Z{i}=Z{i}-diag(diag(Z{i}));
        % Z{i}=max(Z{i},0);
    end

    % TODO: for debug
    % print_matrix(J,S,Z,F,svp,mu);

    % update E_i
    Ek=E;
    [E]=cellfun(@updateE,X,S,Y,cmu,clambda,'UniformOutput',false);
    
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
        fprintf(1,'iter %d,mu is %e,ResidualX is %e,changeZJ is %e,changeZS is %e\n',iter,mu,changeX,changeZJ,changeZS);
        for i=1:k
            fprintf(1,'svp%d %d,',i,svp{i});
        end
        fprintf(1,'\n');
    end
    if changeX <= epsilon && changeZJ <= epsilon && changeZS <= epsilon
    % if changeX <= epsilon && gap <=epsilon2 && changeZJ <= epsilon && changeZS <= epsilon
        convergenced=true;
        fprintf(2,'convergenced, iter is %d\n',iter);
        fprintf(2,'iter %d,mu is %f,ResidualX is %e,changeZJ is %e,changeZS is %e\n',iter,mu,changeX,changeZJ,changeZS);
        for i=1:k
            fprintf(1,'svp%d %d,',i,svp{i});
        end
        fprintf(1,'\n');
        save_matrix(J,S,Z,F);
    end
    % update multipliers
    [Y]=cellfun(@updateY,Y,cmu,Xc,'UniformOutput',false);
    [W]=cellfun(@updateW,W,cmu,ZJc,'UniformOutput',false);
    [V]=cellfun(@updateV,V,cmu,ZSc,'UniformOutput',false);
    % update parameters
    % if gap < epsilon2
        mu=min(rho*mu,max_mu);
    % end
    iter=iter+1;
end
toc


% Jk{i}=J{i};
% [JT,svpt,svt]=singular_value_shrinkage_acc(Z{i}+W{i}/mu,1/mu,sv{i});
% J{i}=JT;
% svp{i}=svpt;
% sv{i}=svt;
function [J, svp, sv] = updateJ(Z,W,mu,sv)
    % [J,svp,sv]=singular_value_shrinkage_acc(Z+W/mu,1/mu,sv);
    [J,svp]=singular_value_shrinkage(Z+W/mu,1/mu); % TODO: sometimes PROPACK is slower than full svd, and sometimes it will throw the following error
    % J = (J - min(J(:))) ./ (max(J(:))-min(J(:))); % TODO: 有可能除0了
    % J=J-diag(diag(J));
    % J=max(J,0);
    % Error using vertcat
    % CAT arguments dimensions are not consistent.

    % Error in lansvd (line 228)
    % [S,bot] = bdsqr(diag(B),[diag(B,-1); resnrm]);

    % Error in singular_value_shrinkage_acc (line 5)
    % [U S V] = lansvd(X, sv, 'L');


% S{i}=invx{i}*(xtx{i}-X{i}'*E{i}+Z{i}+(X{i}'*Y{i}+V{i}-W{i})/mu);
function [S] = updateS(invx,xtx,X,E,Z,Y,V,W,mu)
    S=invx*(xtx-X'*E+Z+(X'*Y+V-W)/mu);
    % S = (S - min(S(:))) ./ (max(S(:))-min(S(:)));
    % S=S-diag(diag(S));
    % S=max(S,0);

% F{i}=(J{i}+S{i}-(W{i}+V{i})*mu)/2;
function [F] = updateF(J,S,W,V,mu)
    % add normalize
    T1=J-W/mu;
    T2=S-V/mu;
    % TODO: normalize not work
    % T1=mnormalize_col(T1);
    % T2=mnormalize_col(T2);
    % F=(J+S-(W+V)*mu)/2; % TODO: fix bug, not *mu, should be /mu
    F=(T1+T2)/2;
    % F = (F - min(F(:))) ./ (max(F(:))-min(F(:)));

% M{i}=reshape(F{i}',1,n*n);
function [M] = updateM(F)
    n=length(F);
    M=reshape(F',1,n*n);

% E{i}=l21(X{i}-X{i}*S{i}+Y{i}/mu,lambda/(2*mu));  % bug fixed, parameter should be lambda/(2*mu), not lambda/mu
function [E] = updateE(X,S,Y,mu,lambda)
    % E=l21(X-X*S+Y/mu,lambda/(2*mu)); % TODO: why 2*mu, I believe it should be mu!
    E=l21(X-X*S+Y/mu,lambda/(mu));

% Xc{i}=X{i}-X{i}*S{i}-E{i}; 
% ZJc{i}=Z{i}-J{i};
% ZSc{i}=Z{i}-S{i};
% vals(i)=norm(Xc{i},'fro')/Xf{i};
% vals(i)=norm(ZJc{i},'fro')/Xf{i};
% vals(i)=norm(ZSc{i},'fro')/Xf{i};

% vals(i)=norm(Zk{i}-Z{i},'fro')/Xf{i};
% vals(i)=norm(Jk{i}-J{i},'fro')/Xf{i};
% vals(i)=norm(Sk{i}-S{i},'fro')/Xf{i};
% vals(i)=norm(Ek{i}-E{i},'fro')/Xf{i};
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


% Y{i}=Y{i}+mu*Xc{i};
% W{i}=W{i}+mu*ZJc{i};
% V{i}=V{i}+mu*ZSc{i};

function [Y] = updateY(Y,mu,Xc)
    Y=Y+mu*Xc;
    
function [W] = updateW(W,mu,ZJc)
    W=W+mu*ZJc;

function [V] = updateV(V,mu,ZSc)
    V=V+mu*ZSc;

function result = errorfun(S, varargin)
   warning(S.identifier, S.message);
   result = NaN;

function [] = save_matrix(J,S,Z,F)
    k=length(J);
    % TODO: for debug
    for i=1:k
        h=figure('Visible', 'off');
        imagesc(S{i});
        colormap(gray);
        axis equal;
        saveas(h,['S' num2str(i) '.png']);
    end
    for i=1:k
        h=figure('Visible', 'off');
        imagesc(Z{i});
        colormap(gray);
        axis equal;
        saveas(h,['Z' num2str(i) '.png']);
    end
    for i=1:k
        h=figure('Visible', 'off');
        imagesc(J{i});
        colormap(gray);
        axis equal;
        saveas(h,['J' num2str(i) '.png']);
    end
    for i=1:k
        h=figure('Visible', 'off');
        imagesc(F{i});
        colormap(gray);
        axis equal;
        saveas(h,['F' num2str(i) '.png']);
    end

function [] = print_matrix(J,S,Z,F,svp,mu)
    % TODO: for debug
    k=length(J);
    for i=1:k
        [aa bb1 cc]=svd(Z{i});
        [aa bb2 cc]=svd(S{i});
        [aa bb3 cc]=svd(F{i});
        fprintf(1,'Z%d rank %d,S%d rank %d,F%d rank %d,J%d rank is %d\n',i,length(find(diag(bb1)>1/mu)),i,length(find(diag(bb2)>1/mu)),i,length(find(diag(bb3)>1/mu)),i,svp{i});
    end
