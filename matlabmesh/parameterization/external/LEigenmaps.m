% --- leigs function for Laplacian eigenmap.
% Written by Belkin & Niyogi, 2002.
function [Y] = LEigenmaps(DATA, K, d) 
n = size(DATA,1);
A = sparse(n,n);
step = 100;  
for i1=1:step:n    
    i2 = i1+step-1;
    if (i2> n) 
      i2=n;
    end;
    XX= DATA(i1:i2,:);  
    dt = L2_distance(XX',DATA',0);
    [Z,I] = sort ( dt,2);
    for i=i1:i2
      for j=2:K+1
	        A(i,I(i-i1+1,j))= Z(i-i1+1,j); 
	        A(I(i-i1+1,j),i)= Z(i-i1+1,j); 
      end;    
    end;
end;
W = A;
[A_i, A_j, A_v] = find(A);  % disassemble the sparse matrix
for i = 1: size(A_i)  
    W(A_i(i), A_j(i)) = 1;
end;
D = sum(W(:,:),2);   
L = spdiags(D,0,speye(size(W,1)))-W;
opts.tol = 1e-9;
opts.issym=1; 
opts.disp = 0; 
[E,V] = eigs(L,d+1,'sm',opts);

Y = E(:,1:d);