
mesh = readMesh('patch.obj');
mesh = readMesh('patch2.obj', 'n');

% 1) choose boundary embedding
boundaryUV = embedBoundary( mesh, 'circle' );

% 2) compute weights and embed

weights = makeOneRingWeights(mesh, 'invdist');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

weights = makeOneRingWeights(mesh, 'uniform');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

weights = makeOneRingWeights(mesh, 'gaussian');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');


weights = makeOneRingWeights(mesh, 'meanvalue');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

weights = makeOneRingWeights(mesh, 'dcp');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

% save parameterized mesh
writeMesh(mesh,'temp.obj');