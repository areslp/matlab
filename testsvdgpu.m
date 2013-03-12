N=4096;
A=rand(N);
tic;
[u,s,v]=culasvd(A);
toc
tic
[u,s,v]=svd(A);
toc