function [ deformed_mesh ] = deformLaplacian( mesh, vertTrans, consP )
%[ deformed_mesh ] = deformLaplacian( mesh, vertTrans, consP )
%
%   vertTrans:     per-vertex transformation matrices
%                   - Vx9 matrix with rows  [x11,x12,x12,x21,...,x33]
%   consP:          rows of [vtx_i, x, y, z, weight_i]

N = numel(mesh.vidx);

% construct Laplacian differential operator  ("Ls" in [BS07] survey)
Ls = makeOneRingWeights(mesh, 'cotan', 0) / 2;
Ls(1:N+1:N*N) = -sum(Ls,2); 

% construct per-vertex weight/mass matrix ('M' in [BS07] survey)
areas = vertexArea(mesh, [], 'mixed');
Minv = sparse(mesh.vidx,mesh.vidx,  areas.^-1,  N,N);

% Laplacian operator with weights taken into account
L = Minv*Ls;

% construct rotated laplacians
delta = L*mesh.v;
Rdelta = zeros(size(mesh.v));
for i = 1:N
    trans = reshape( vertTrans(i,:), 3,3 );
    l = delta(i,:);
    Rdelta(i,:) = mvmul(trans,l);
end

% construct normal equations
Msys = Ls * Minv * Ls;
RHS = Ls * Rdelta;

% add soft constraints
consi = consP(:,1);
consv = consP(:,2:4);
consw = consP(:,5);
[Mcons,RHScons] = softConstrain(Msys,RHS,consi,consv,consw);

% solve and unpack
X = Mcons \ RHScons;
deformed_mesh = mesh;
deformed_mesh.v = [X(:,1), X(:,2), X(:,3)];
deformed_mesh.n = estimateNormal(deformed_mesh);

end
