function [ Mout ] = vdot( M, M2 )
%VDOT returns vector of dot products of rows of M and M2. If M2 is a single
% row, expands to # of rows in M1

if size(M2,1) == 1
    n = size(M,1);
    Mout = sum( M .* repmat(M2, size(M,1), 1), 2 );
else
    Mout = sum( M .* M2, 2 );
end
