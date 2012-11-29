mesh = readMesh('patch2.obj');
mesh.n = estimateNormal(mesh, mesh.vidx);
plotMesh(mesh);

faired = fair_Laplacian(mesh);
plotMesh(faired);



mesh = readMesh('fill.obj');
plotMesh(mesh,'nf');
mesh.n = estimateNormal(mesh, mesh.vidx);

faired = fair_LBH(mesh, 10, 0.0001, 1);
plotMesh(faired);

faired2 = fair_LBH(faired2, 250, 0.00001);


faired = mesh;
faired.n = estimateNormal(faired, faired.vidx);

faired = fair_LBH(faired, 10);
faired.n = estimateNormal(faired, faired.vidx);
plotMesh(faired);




