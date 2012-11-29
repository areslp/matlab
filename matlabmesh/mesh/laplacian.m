function [ L ] = laplacian( mesh, vertices, mode )
%MEANCURV estimate vertex laplacians
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   - If vertices empty or undefined compute laplacians for entire mesh 
%   - Valid modes are:
%     'uniform' - uniform weights
%     'cotan' - cotangent weights


if ~ exist('mode', 'var')
    mode = 'uniform';
end

if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = 1:size(mesh.v,1);
end

size(vertices)
mode

n = numel(vertices);
dim = size(mesh.v,2);
   
L = zeros(n,dim);
for i = 1:n
    if strcmp(mode, 'uniform')
        L(i,:) = laplacian_uniform( mesh, vertices(i) );
    elseif strcmp(mode, 'cotan')
        L(i,:) = laplacian_cotan( mesh, vertices(i) );
    end
end

end


function [ L ] = laplacian_cotan( mesh, v )
    [ot,ov] = onering(mesh, v);
    w = cotanWeights(mesh, v);
    wj = nonzeros(w(1,ov));
    wj = wj / sum(wj);
    qj = mesh.v(ov,:);
    c = sum( vmul(qj, wj), 1);
    L = c - mesh.v(v,:);
end

function [ L ] = laplacian_uniform( mesh, v )
    [ot,ov] = onering(mesh, v);
    qj = mesh.v(ov,:);
    c = sum(qj,1) / numel(ov);
    L = c - mesh.v(v,:);
end