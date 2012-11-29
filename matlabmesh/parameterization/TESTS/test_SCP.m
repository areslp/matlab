%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick an input mesh to read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% portion of a sphere
mesh = readMesh('patch.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

mesh = readMesh('patch2.obj', 'n');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');


%% parameterize

mesh.uv = embedSCP( mesh, 'fiedler' );
plotMesh(mesh, 'uefb');

mesh.uv = embedSCP( mesh, 'generalized' );
plotMesh(mesh, 'uefb');

mesh.uv = embedSCP( mesh, 'robust' );
plotMesh(mesh, 'uefb');

mesh.uv = embedSCP( mesh, 'lle' );
plotMesh(mesh, 'uefb');




