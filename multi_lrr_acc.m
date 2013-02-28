function [ZZ,Z,E] = multi_lrr(X,lambda,alpha)
% implement the algorithm described in paper "Multi-task Low-rank Affinity Pursuit for Image Segmentation"
% for test use, hard code the affinity matrix number
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
epsilon=1e-4;
epsilon2=1e-5; % must be small!

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

MAX_ITER=1000;
iter=0;
convergenced=false;



tic
while ~convergenced
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        break;
    end
    % tic
    % update J_i
    for i=1:k
        Jk{i}=J{i};
        [JT,svpt,svt]=singular_value_shrinkage_acc(Z{i}+W{i}/mu,1/mu,sv{i});
        % [JT,svpt]=singular_value_shrinkage(Z{i}+W{i}/mu,1/mu);
        J{i}=JT;
        svp{i}=svpt;
        sv{i}=svt;
    end
    % t=toc;
    % fprintf(1,'singular_value_shrinkage takes %f\n',t);
    % update S_i
    % tic
    % for i=1:k
        % Sk{i}=S{i};
        % S{i}=invx{i}*(xtx{i}-X{i}'*E{i}+Z{i}+(X{i}'*Y{i}+V{i}-W{i})/mu);
    % end
    Sk=S;
    cmu=cell(k,1);
    cmu(1:k)={mu};
    S=cellfun(@updateS,invx,xtx,X,E,Z,Y,V,W,cmu,'UniformOutput',false,'ErrorHandler',@errorfun);
    % t=toc;
    % fprintf(1,'update S takes %f\n',t);
    % update Z
    % tic;
    for i=1:k
        Zk{i}=Z{i};
        F{i}=(J{i}+S{i}-(W{i}+V{i})*mu)/2;
    end
    % t=toc;
    % fprintf(1,'update Z takes %f\n',t);
    % update ZZ
    % tic;
    M=cell(k,1);
    for i=1:k
        M{i}=reshape(F{i}',1,n*n);
    end
    MM=zeros(k,n*n);
    for i=1:k
        MM(i,:)=M{i};
    end
    ZZ=l21(MM,alpha/(2*mu));
    % update Z_i
    for i=1:k
        Zk{i}=Z{i};
        Z{i}=reshape(ZZ(i,:),n,n)';
    end
    % t=toc;
    % fprintf(1,'update ZZ and Z_i takes %f\n',t);
    % update E_i
    % tic;
    for i=1:k
        Ek{i}=E{i};
        E{i}=l21(X{i}-X{i}*S{i}+Y{i}/mu,lambda/(2*mu));  % bug fixed, parameter should be lambda/(2*mu), not lambda/mu
    end
    % t=toc;
    % fprintf(1,'update E takes %f\n',t);
    
    % tic;
    for i=1:k
       Xc{i}=X{i}-X{i}*S{i}-E{i}; 
       ZJc{i}=Z{i}-J{i};
       ZSc{i}=Z{i}-S{i};
    end
    % check convergence
    % find the max error in multi views
    vals=zeros(k,1);
    for i=1:k
        vals(i)=norm(Xc{i},'fro')/Xf{i};
    end
    changeX=max(vals);
    for i=1:k
        vals(i)=norm(ZJc{i},'fro')/Xf{i};
    end
    changeZJ=max(vals);
    for i=1:k
        vals(i)=norm(ZSc{i},'fro')/Xf{i};
    end
    changeZS=max(vals);
    for i=1:k
        vals(i)=norm(Zk{i}-Z{i},'fro')/Xf{i};
    end
    changeZ=max(vals); 
    for i=1:k
        vals(i)=norm(Jk{i}-J{i},'fro')/Xf{i};
    end
    changeJ=max(vals);
    for i=1:k
        vals(i)=norm(Sk{i}-S{i},'fro')/Xf{i};
    end
    changeS=max(vals);
    for i=1:k
        vals(i)=norm(Ek{i}-E{i},'fro')/Xf{i};
    end
    changeE=max(vals);
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
    for i=1:k
        Y{i}=Y{i}+mu*Xc{i};
    end
    for i=1:k
        W{i}=W{i}+mu*ZJc{i};
    end
    for i=1:k
        V{i}=V{i}+mu*ZSc{i};
    end
    % update parameters
    if gap < epsilon2
        mu=min(rho*mu,max_mu);
    end
    iter=iter+1;
    % t=toc;
    % fprintf(1,'other things takes %f\n',t);
end
toc
% S{i}=invx{i}*(xtx{i}-X{i}'*E{i}+Z{i}+(X{i}'*Y{i}+V{i}-W{i})/mu);
function [S] = updateS(invx,xtx,X,E,Z,Y,V,W,mu)
    S=invx*(xtx-X'*E+Z+(X'*Y+V-W)/mu);

function result = errorfun(S, varargin)
   warning(S.identifier, S.message);
   result = NaN;

