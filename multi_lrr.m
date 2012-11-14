function [Z,Z1,Z2,Z3,E1,E2,E3] = multi_lrr(X1,X2,X3,lambda,alpha)
% implement the algorithm described in paper "Multi-task Low-rank Affinity Pursuit for Image Segmentation"
% for test use, hard code the affinity matrix number
[m,n]=size(X1);
% initial matrix
E1=zeros(m,n);
E2=zeros(m,n);
E3=zeros(m,n);

J1=zeros(n,n);
J2=zeros(n,n);
J3=zeros(n,n);

S1=zeros(n,n);
S2=zeros(n,n);
S3=zeros(n,n);

Z1=zeros(n,n);
Z2=zeros(n,n);
Z3=zeros(n,n);

W1=zeros(n,n);
W2=zeros(n,n);
W3=zeros(n,n);

Y1=zeros(m,n);
Y2=zeros(m,n);
Y3=zeros(m,n);

V1=zeros(n,n);
V2=zeros(n,n);
V3=zeros(n,n);

Z=zeros(3,n*n);

% parameters
mu=1e-6;
max_mu=10^10;
rho=1.1;
epsilon=1e-6;
epsilon2=1e-2;

% pre caculate matrix value
xtx1=X1'*X1;
inv_x1=inv(xtx1+eye(n));
xtx2=X2'*X2;
inv_x2=inv(xtx2+eye(n));
xtx3=X3'*X3;
inv_x3=inv(xtx3+eye(n));
X1f=norm(X1,'fro');
X2f=norm(X2,'fro');
X3f=norm(X3,'fro');

MAX_ITER=1000;
iter=0;
convergenced=false;

while ~convergenced
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        break;
    end
    % update J_i
    J1b=J1;
    J2b=J2;
    J3b=J3;
    J1=singular_value_shrinkage(Z1+W1/mu,1/mu);
    J2=singular_value_shrinkage(Z2+W2/mu,1/mu);
    J3=singular_value_shrinkage(Z3+W3/mu,1/mu);
    % update S_i
    S1b=S1;
    S2b=S2;
    S3b=S3;
    S1=inv_x1*(xtx1-X1'*E1+Z1+(X1'*Y1+V1-W1)/mu);
    S2=inv_x2*(xtx2-X2'*E2+Z2+(X2'*Y2+V2-W2)/mu);
    S3=inv_x3*(xtx3-X3'*E3+Z3+(X3'*Y3+V3-W3)/mu);
    % update Z
    Zb=Z;
    F1=(J1+S1-(W1+V1)*mu)/2;
    F2=(J2+S2-(W2+V2)*mu)/2;
    F3=(J3+S3-(W3+V3)*mu)/2;
    M=[reshape(F1',1,n*n);reshape(F2',1,n*n);reshape(F3',1,n*n)];
    Z=l21(M,alpha/(2*mu));
    % update Z_i
    Z1b=Z1;
    Z2b=Z2;
    Z3b=Z3;
    Z1=reshape(M(1,:),n,n)';
    Z2=reshape(M(2,:),n,n)';
    Z3=reshape(M(3,:),n,n)';
    % update E_i
    E1b=E1;
    E2b=E2;
    E3b=E3;
    E1=l21(X1-X1*S1+Y1/mu,lambda/mu); 
    E2=l21(X2-X2*S2+Y2/mu,lambda/mu); 
    E3=l21(X3-X3*S3+Y3/mu,lambda/mu); 

    % temporary matrix
    X1_c=X1-X1*S1-E1;
    X2_c=X2-X2*S2-E2;
    X3_c=X3-X3*S3-E3;

    ZJ1_c=Z1-J1;
    ZJ2_c=Z2-J2;
    ZJ3_c=Z3-J3;

    ZS1_c=Z1-S1;
    ZS2_c=Z2-S2;
    ZS3_c=Z3-S3;

    % check convergence
%     changeX=max(max(norm(X1_c,'fro')/X1f,norm(X2_c,'fro')/X2f),norm(X3_c,'fro')/X3f);
%     changeZJ=max(max(norm(ZJ1_c,'fro')/X1f,norm(ZJ2_c,'fro')/X2f),norm(ZJ3_c,'fro')/X3f);
%     changeZS=max(max(norm(ZS1_c,'fro')/X1f,norm(ZS2_c,'fro')/X2f),norm(ZS3_c,'fro')/X3f);
    
    changeX=max(max(norm(X1_c,'fro')/X1f,norm(X2_c,'fro')/X2f),norm(X3_c,'fro')/X3f);
    changeZJ=max(max(norm(ZJ1_c,'fro')/X1f,norm(ZJ2_c,'fro')/X2f),norm(ZJ3_c,'fro')/X3f);
    changeZS=max(max(norm(ZS1_c,'fro')/X1f,norm(ZS2_c,'fro')/X2f),norm(ZS3_c,'fro')/X3f);
    fprintf(1,'mu is %f,changeX is %f,changeZJ is %f,changeZS is %f\n',mu,changeX,changeZJ,changeZS);
    if changeX <= epsilon && changeZJ <= epsilon && changeZS <= epsilon
        convergenced=true;
        fprintf(1,'convergenced, iter is %d\n',iter);
    end
    % if changeX <= epsilon
        % convergenced=true;
        % fprintf(1,'convergenced, iter is %d\n',iter);
    % end

    % update multipliers
    Y1=Y1+mu*(X1_c);
    Y2=Y2+mu*(X2_c);
    Y3=Y3+mu*(X3_c);

    W1=W1+mu*(ZJ1_c);
    W2=W2+mu*(ZJ2_c);
    W3=W3+mu*(ZJ3_c);

    V1=V1+mu*(ZS1_c);
    V2=V2+mu*(ZS2_c);
    V3=V3+mu*(ZS3_c);
    % update parameters
    changeZ=max(max(norm(Z1-Z1b,'fro')/X1f,norm(Z2-Z2b,'fro')/X2f),norm(Z3-Z3b,'fro')/X3f);
    changeJ=max(max(norm(J1-J1b,'fro')/X1f,norm(J2-J2b,'fro')/X2f),norm(J3-J3b,'fro')/X3f);
    changeS=max(max(norm(S1-S1b,'fro')/X1f,norm(S2-S2b,'fro')/X2f),norm(S3-S3b,'fro')/X3f);
    changeE=max(max(norm(E1-E1b,'fro')/X1f,norm(E2-E2b,'fro')/X2f),norm(E3-E3b,'fro')/X3f);

    tmp=[changeZ changeJ changeS changeE];
    fprintf(1,'gap between two iteration is %f\n',max(tmp));
    if max(tmp) < epsilon2
        mu=min(rho*mu,max_mu);
    end
    
    iter=iter+1;
    fprintf(1,'iter is %d\n',iter);
end
