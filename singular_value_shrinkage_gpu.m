function [ Y,svp ] = singular_value_shrinkage_gpu( X, tau )
[U,S,V]=culasvd(X);
% [U,S,V]=lansvd(X);
S = diag(S);
svp = length(find(S > tau));
% fprintf(1,'svp:%d\n',svp);
S = S(1:svp)-tau;
U = U(:, 1:svp);
V = V(:, 1:svp);
Y = U*diag(S)*V';
end


