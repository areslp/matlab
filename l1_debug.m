function [result] = l1_debug(xi,X)
% return the dot value of xi and the columns of X
[m n]=size(X);
epsilon=1e-4;
% norm(xi)
assert(abs(norm(xi)-1)<epsilon); %xi and the columns of X must be normalized
result=zeros(1,n);
for i=1:n
    cX=X(:,i);
    % norm(cX)
    assert(abs(norm(cX)-1)<epsilon);
    result(i)=xi'*cX;
end
end
