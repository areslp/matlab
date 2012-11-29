% pick an input mesh to read

mesh = readMesh('patch.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

mesh = readMesh('patch2.obj', 'n');
plotMesh(mesh,'efb');

mesh = readMesh('bunny.obj');
plotMesh(mesh,'efb');


% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
plotMesh(mesh,'efb');

mesh = readMesh('doghead.obj');
plotMesh(mesh,'efb');

mesh = readMesh('doghead_cut2.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

mesh = readMesh('4bump.obj');
plotMesh(mesh,'efb');

mesh = readMesh('doghead_basehole.obj');
plotMesh(mesh,'efb');



% compute DNCP (discrete natural-boundary conformal) parameterization
mesh.u = embedDNCP(mesh);
figure;plotMesh(mesh, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
fprintf('[DNCP ] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  


% compute SCP (spectral conformal) parameterization
mesh.u = embedSCP(mesh, 'fiedler');
plotMesh(mesh, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
fprintf('[SCP Fiedler] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  

% compute SCP (spectral conformal) parameterization
mesh.u = embedSCP(mesh, 'generalized');
plotMesh(mesh, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
fprintf('[SCP Generalized] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  




tmp = mesh;
tmp.v = embedSCP(mesh, 'fiedler');
plotMesh(tmp, 'efb');


% compare uniform/DCP/DAP parameterizations w/ circle boundary
write_meshes = 0;   % set to 1 to write parameterized meshes

boundaryUV = embedBoundary( mesh, 'circle' );

weights = makeOneRingWeights(mesh, 'uniform');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
fprintf('[Uniform  ] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  
if write_meshes
    writeMesh(mesh, '_param_uniform.obj');
end

weights = makeOneRingWeights(mesh, 'dcp');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
fprintf('[Conformal] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  
if write_meshes
    writeMesh(mesh, '_param_conformal.obj');
end

weights = makeOneRingWeights(mesh, 'dap');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
fprintf('[Authalic ] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  
if write_meshes
    writeMesh(mesh, '_param_authalic.obj');
end



% check that planar mesh is reproduced by DCP/DAP

tmp = mesh;
tmp.v = [mesh.u(:,1), mesh.u(:,2), zeros(numel(mesh.vidx),1)];
boundaryUV = embedBoundary( tmp, 'copy' );
plotMesh(tmp,'efb');
weights = makeOneRingWeights(tmp, 'dcp');
tmp.u = embedInterior(tmp, boundaryUV, weights);
plotMesh(mesh, 'uefb');
plotMesh(tmp, 'uefb');
[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(tmp, 1);
fprintf('[Conformal] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  
Edcp = sum(sum(mesh.u-tmp.u));
weights = makeOneRingWeights(tmp, 'dap');
tmp.u = embedInterior(tmp, boundaryUV, weights);
plotMesh(mesh, 'uefb');
Edap = sum(sum(mesh.u-tmp.u));
fprintf(' DCP error: %.9f   DAP error: %.9f\n', Edcp, Edap);


