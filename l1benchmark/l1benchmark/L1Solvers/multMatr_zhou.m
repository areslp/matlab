function y=multMatr_zhou(A,x)
%y=A*x
%
% Call: y=multMatr(A,x)

k = size(A,2);
y=A*x(1:k) + x(k+1:end);