function y=multMatrAdj_zhou(A,x)
%y=A'*x
%
% Call: y=multMatrAdj(A,x)

y=([x'*A x'])';