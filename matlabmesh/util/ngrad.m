function [ df ] = ngrad( f, x, dx )
%NGRAD numerical gradient approximation
%   RMSTODO: check if this actually works, let dx be a vector
n = size(x,1);
val = f(x);
for i = 1:n
    xi = x;
    xi(i) = xi(i) + dx;
    df(i) = val - f(xi);
end
