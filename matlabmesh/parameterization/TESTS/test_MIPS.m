% read a test mesh
mesh = readMesh('patch.obj');
mesh = readMesh('patch2.obj', 'n');

plotMesh(mesh, 'vefb');


% compute initial embedding on circle using uniform weights
boundaryUV = embedBoundary( mesh, 'circle' );
weights = makeOneRingWeights(mesh, 'uniform');
uvmesh = mesh;
uvmesh.u = embedInterior(mesh, boundaryUV, weights);

plotMesh(uvmesh, 'uefb');


% minimize MIPS energy using K optimization steps
K = 2000;
uvmesh = embedMIPS( mesh, K );
plotMesh(uvmesh, 'uefb');

