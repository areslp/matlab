function [Z] = l1_graph( X, lambda)

features=X;

for i=1:size(features,2)
   features(:,i)= features(:,i)/norm(features(:,i));
end

f_dim=size(features,1);
num=size(features,2);

f=fopen('spmatrix.txt','w');
for i=1:num
    X=[features(:,1:i-1) features(:,i+1:num)];
    B=[X];
    xi=features(:,i);

    % alpha=LassoNonNegativeSquared(B,xi,lambda);
    % [alpha, nIter, timeSteps, errorSteps] = SolveDALM(B, xi, 'lambda',lambda, 'stoppingcriterion',1);
    opt.rho=lambda;
    opt.nonneg=1;
    opt.tol=1e-3;
    alpha=yall1(B,xi,opt);
    
    %remove the non-zero alpha on I
    assert(length(alpha)==num-1);
    alpha=alpha(1:num-1);
    % if length(find(alpha~=0))==0
        % fprintf(1,'alpha is 0, i is %d\n',i);
        % displayPatches(B);
        % dot_result=l1_debug(xi,B)
        % pause;
    % end

    [r,c,v]=find(alpha);
    
    for ii=1:length(r)
        row=i;
        col=r(ii);
        if col<i
        else
            col=col+1;
        end
        % col=c(ii);  all 1
        val=abs(v(ii));
        fprintf(f,'%d %d %f\n',row,col,val);
    end
end
fclose(f);

load spmatrix.txt
A=spconvert(spmatrix);
num=length(A);
if size(A,1)~=size(A,2)
    A(num,num)=0;
end

Z=A;

