
mesh = readMesh('patch.obj');
mesh = readMesh('patch2.obj', 'n');
mesh = readMesh('doghead.obj');
mesh = readMesh('bump_ref.obj');

%plotMesh(mesh,'vefbn');
plotMesh(mesh,'efb');

% 1) choose boundary embedding

boundaryUV = embedBoundary( mesh, 'circle' );


% 2) compute weights

% optimal3 weights result in foldovers, etc, if
% nbrhoods are too large (negative weights?)
tic;
options.weightmode = 'optimal3';
options.nbrtype = 'k';
options.nbrsize = 8;
weights = makeExpMapWeights(mesh.v, mesh.n, options);
toc;

% optimal2 weights can be negative, but don't have above problem...
tic;
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 8;
weights = makeExpMapWeights(mesh.v, mesh.n, options);
toc;

% slightly better with geodesic-ball neighbourhoods...
tic;
edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);
toc;


% 3) embed interior

mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');



% save parameterized mesh
writeMesh(mesh,'temp.obj');






% embed using LLE

options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 10;
options.regtol = 1e-5;
weights = makeExpMapWeights(mesh.v, mesh.n, options);

[mesh.u,xyz] = embedLLE(mesh.v, weights);
plotMesh(mesh, 'uefb');

tmpmesh = mesh;
tmpmesh.v = xyz;
plotMesh(tmpmesh, 'efb');


mesh.u = lle(mesh.v', 5, 2)';
plotMesh(mesh, 'uefb');


