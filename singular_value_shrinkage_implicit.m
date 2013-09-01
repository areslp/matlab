<<<<<<< HEAD
function [ U,S,V ] = singular_value_shrinkage( X, tau )
% tau
=======
function [ U,S,V,svp ] = singular_value_shrinkage( X, tau )
>>>>>>> 4e2789888db39b134173590ff8bfd29ee6c32a79
%SINGULAR_VALUE_SHRINKAGE Summary of this function goes here
%   Detailed explanation goes here
% tic;
[U,S,V]=svd(X,'econ');
% [U,S,V]=svd(X);
% toc
% [U,S,V]=lansvd(X);
S = diag(S);
svp = length(find(S > tau));
% fprintf(1,'svp:%d\n',svp);
S = S(1:svp)-tau;
% S=diag(wthresh(diag(S),'s',tau)); % 在实现里它只考虑>tau的情况
U = U(:, 1:svp);
V = V(:, 1:svp);
% Y = U*diag(S)*V';
% t=toc;
% fprintf(1,'svd shrinkage takes:%f\n',t);
end


