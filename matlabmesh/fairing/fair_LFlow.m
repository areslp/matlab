function [ mesh ] = fair_LFlow( init_mesh, speed, iters )
%FAIR_LFLOW try to fair mesh using simple laplacian flow

mesh = init_mesh;

n = size(mesh.v,1);
Bi = mesh.vidx( mesh.isboundaryv == 1 );  % boundary verts 
Ii = mesh.vidx( mesh.isboundaryv == 0 );  % interior verts 

if ~exist('speed','var')
    speed = 0.1;
end
if ~exist('iters','var')
    iters = 10;
end

for iter = 1:iters
    fprintf('iter %d\n', iter);

    L = laplacian(mesh,Ii,'uniform');
    mesh.v(Ii,:) = mesh.v(Ii,:) - speed * L;

    plotMesh(mesh,'efbn');
    drawnow;
end

end
