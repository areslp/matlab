function [ uv ] = embedInterior( mesh, boundaryUV, weights, constraints )
% [ uv ] = embedInterior( mesh, boundaryUV, weights, constraints )
% compute 2D embedding (parameterization) of mesh
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
% embed interior verts of mesh in 2D based on boundary
%  UV embedding and interior neighbourhood weights
% RMSTODO: solving bigger matrix than necessary. Also
%   solving it twice. Should be using LU decomposion...

n = size(mesh.v,1);
uv = zeros(n,2);
uv( boundaryUV(:,1), :) = [boundaryUV(:,2), boundaryUV(:,3)];

M = spalloc(n, n, nnz(weights) + n);
rhsu = zeros(n,1);  rhsv = zeros(n,1);
for v = 1:n
    [i,j, w] = find(weights(v,:));

    bIsBoundary = ~isempty( find(mesh.loops{1} == v) );
    if bIsBoundary
        M(v,v) = 1;
        rhsu(v) = uv(v,1);
        rhsv(v) = uv(v,2);
    else
%         if min(w) < 0
%             fprintf('min weight < 0 for vtx %d\n',v);
%         end
        M(v, j) = -w;
        M(v, v) = sum(w);
    end
end


% apply linear position constraints by clearing rows
if exist('constraints','var')
    for ci = 1:size(constraints,1)
        v = constraints(ci,1);
        M(v,:) = 0;
        M(v,v) = 1;
        rhsu(v) = constraints(ci,2);
        rhsv(v) = constraints(ci,3);
    end
end

solveu = M \ rhsu;
solvev = M \ rhsv;

uv = [solveu, solvev];

