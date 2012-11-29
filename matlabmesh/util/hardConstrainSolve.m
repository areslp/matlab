function [ X ] = hardConstrainSolve( M, B, consi, consv )
% [ X ] = hardConstrainSolve( M, RHS, consi, consv )
%   solve linear system M.X - B with hard constraints  (automatically handles rewriting indices)
%     consi = row indices of values to constrain
%     consv = constrained values

[nrows,ncols] = size(M);
keepi = setdiff(1:nrows,consi);

[Mc,Bc] = hardConstrain(M,B,consi,consv, 1);
Xc = Mc \ Bc;

X = zeros(size(B));
X(consi,:) = consv;
X(keepi,:) = Xc;

