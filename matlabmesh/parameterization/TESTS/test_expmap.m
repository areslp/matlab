% read a test mesh
mesh = readMesh('patch.obj');
nVtx = 63;   % this is the vertex that expmap will be centered around

mesh = readMesh('patch2.obj', 'n');
nVtx = 324;

plotMesh(mesh, 'vefbn');


% run expmap

uvmesh = mesh;

[uvmesh.u,geo_dists] = expmap(mesh.v, mesh.n, nVtx);

plotMesh(uvmesh, 'uefb');



% run expmap with K-nbr stop criteria
uvmesh = mesh;
options = [];
options.stopcrit = 'k';
options.stopparam = 60;
[uvmesh.u,geo_dists] = expmap(mesh.v, mesh.n, nVtx, options);
plotMesh(uvmesh, 'uefb');


% run expmap with geodesic-ball stop criteria
edgelen_avg = edgeStats(mesh.v, mesh.e);
uvmesh = mesh;
options = [];
options.stopcrit = 'maxdist';
options.stopparam = 5*edgelen_avg;
[uvmesh.u,geo_dists] = expmap(mesh.v, mesh.n, nVtx, options);
plotMesh(uvmesh, 'uefb');


% expmap caching test
tic;
[uvmesh.u,geo_dists,options_cache] = expmap(mesh.v, mesh.n, nVtx);
toc;
tic;
[uvmesh.u,geo_dists] = expmap(mesh.v, mesh.n, nVtx, options_cache);
toc;


% this should look like a cone
uvmesh.v = [uvmesh.u(:,1), uvmesh.u(:,2), geo_dists];
plotMesh(uvmesh, 'efb');

% this should look like a smoother cone, except for far-away distortion...
uvmesh.v = [uvmesh.u(:,1), uvmesh.u(:,2), vmag(uvmesh.u)];
plotMesh(uvmesh, 'efb');


% save parameterized mesh
writeMesh(uvmesh,'temp.obj');