%% selection tests

y = @(x,y,z) implicitSphere(x,y,z, [0,0,0.45], 0.5);
y = @(x,y,z) implicitPlane(x,y,z, [0,0,0], [0,0,1]);

fidx = selectMesh(mesh, y, 'f');
subm = subMesh(mesh, fidx);
plotMesh(subm);




%% meshDeform_Lipman04 tests

mesh = readMesh('patch.obj');
handle_vi = 65;
handle_pos = [0.2,0.2,1.0];
roi = mesh.fidx;
[roi_mesh,vmap,fmap] = subMesh(mesh, roi);

mesh = readMesh('hemi_bumpy_lo.obj');
mesh.n = estimateNormal(mesh);
handle_vi = 91;
handle_pos = mesh.v(handle_vi,:) + [2.0,0,1.0];
%handle_pos = mesh.v(handle_vi,:) + [0,0,0.0];
roi = mesh.fidx;
[roi_mesh,vmap,fmap] = subMesh(mesh, roi);

mesh = readMesh('sphere.obj');
handle_vi = 127;
handle_pos = mesh.v(handle_vi,:) + [0.5,0,0.5];
y = @(x,y,z) implicitSphere(x,y,z, mesh.v(handle_vi,:), 0.5);
roi = selectMesh(mesh, y, 'f');
[roi_mesh,vmap,fmap] = subMesh(mesh, roi);
plotMesh(roi_mesh, 'vefb');
handle_vi = vmap( vmap(:,2)==handle_vi , 1 );


% set boundary and handle constraints
wbdry = 100.0;
winterior = 1.0;
Bi = roi_mesh.vidx( roi_mesh.isboundaryv == 1 );  % boundary verts
constraints = [Bi, roi_mesh.v(Bi,:), repmat(wbdry, numel(Bi), 1) ];
constraints = [constraints; handle_vi,  handle_pos(1), handle_pos(2), handle_pos(3), winterior];

% compute deformed roi and update roi vertices
rotation_iters = 0;
deformed_roi = meshDeform_Lipman04(roi_mesh, constraints, rotation_iters);
deformed_mesh = mesh;
deformed_mesh.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);

% plot result, constraints, and constraint results
hold all;
plotMesh(deformed_mesh);
scatter3( constraints(:,2), constraints(:,3), constraints(:,4) );
ci = constraints(:,1);
scatter3( deformed_roi.v(ci,1), deformed_roi.v(ci,2), deformed_roi.v(ci,3) );
hold off;


% compare results with two different parameter settings (shifts second along x)
deformed_roi = meshDeform_Lipman04(roi_mesh, constraints, 0);
deformed_1 = mesh;
deformed_1.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);

deformed_roi = meshDeform_Lipman04(roi_mesh, constraints, 5);
deformed_2 = mesh;
deformed_2.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);
width = mesh.bounds(2,1) - mesh.bounds(1,1);
deformed_2.v = vadd(deformed_2.v,[1.25*width,0,0]);

hold on;
plotMesh(deformed_1);
plotMesh(deformed_2);
hold off;







%% laplacian boundary deformation test


mesh = readMesh('bunny.obj');
mesh.n = estimateNormal(mesh);
bentmesh = readMesh('bunny_bent.obj');
bentmesh.n = estimateNormal(bentmesh);

mesh = readMesh('uniplane.obj');
mesh.n = estimateNormal(mesh);
bentmesh = readMesh('uniplane_bent.obj');
bentmesh.n = estimateNormal(bentmesh);

mesh = readMesh('hemicube.obj');
mesh.n = estimateNormal(mesh);
bentmesh = readMesh('hemicube_bent.obj');
bentmesh.n = estimateNormal(bentmesh);


mesh = readMesh('4bump.obj');
mesh.n = estimateNormal(mesh);
bentmesh = readMesh('4bump_deformed.obj');
bentmesh.n = estimateNormal(bentmesh);


mesh = readMesh('mushroom_lo.obj');
mesh.n = estimateNormal(mesh);
bentmesh = readMesh('mushroom_lo_bent.obj');
bentmesh.n = estimateNormal(bentmesh);



roi = mesh.fidx;
[roi_mesh,vmap,fmap] = subMesh(mesh, roi);
bent_roi = subMesh(bentmesh,roi);

wbdry = 100.0;
Bi = roi_mesh.vidx( roi_mesh.isboundaryv == 1 );  % boundary verts
constraints = [Bi, bent_roi.v(Bi,:), repmat(wbdry, numel(Bi), 1) ];

% compute deformed roi and update roi vertices
rotation_iters = 0;
deformed_roi = meshDeform_Lipman04(roi_mesh, constraints, rotation_iters);
deformed_mesh = mesh;
deformed_mesh.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);

% plot result, constraints, and constraint results
hold all;
plotMesh(deformed_mesh);
scatter3( constraints(:,2), constraints(:,3), constraints(:,4) );
ci = constraints(:,1);
scatter3( deformed_roi.v(ci,1), deformed_roi.v(ci,2), deformed_roi.v(ci,3) );
hold off;


% compare results
orig = mesh;
width = mesh.bounds(2,1) - mesh.bounds(1,1);
origbent = bentmesh;
origbent.v = vadd(origbent.v,[1.25*width,0,0]);


deformed_roi = meshDeform_Lipman04(roi_mesh, constraints, 0);
deformed_1 = mesh;
deformed_1.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);
deformed_1.v = vadd(deformed_1.v,[2.5*width,0,0]);

deformed_roi = meshDeform_Lipman04(roi_mesh, constraints, 3);
deformed_2 = mesh;
deformed_2.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);
deformed_2.v = vadd(deformed_2.v,[3.5*width,0,0]);

deformed_roi = meshDeform_Laplacian(roi_mesh, constraints, bentmesh.n);
deformed_3 = mesh;
deformed_3.v(vmap(:,2),:) = deformed_roi.v(vmap(:,1),:);
deformed_3.v = vadd(deformed_3.v,[4.5*width,0,0]);

hold on;
plotMesh(orig);
plotMesh(origbent);
plotMesh(deformed_1);
plotMesh(deformed_2);
plotMesh(deformed_3);
hold off;




%% poisson mesh reconstruction test 


mesh = readMesh('hemicube.obj');
mesh.n = estimateNormal(mesh);
bentmesh = readMesh('hemicube_bent.obj');
bentmesh.n = estimateNormal(bentmesh);

Bi = mesh.vidx( mesh.isboundaryv == 1 );  % boundary verts

% identity transformation...
deformed_tris = [ mesh.v(mesh.f(:,1),:), mesh.v(mesh.f(:,2),:), mesh.v(mesh.f(:,3),:) ];
boundary_cons = [Bi, mesh.v(Bi,:) ];

tmpmesh = mesh;
tmpmesh.v(43,2) = tmpmesh.v(43,2) + 0.5;
deformed_tris = [ tmpmesh.v(mesh.f(:,1),:), tmpmesh.v(mesh.f(:,2),:), tmpmesh.v(mesh.f(:,3),:) ];
boundary_cons = [Bi, tmpmesh.v(Bi,:) ];

shiftmesh = mesh;
shiftmesh.v = vadd(shiftmesh.v,[0.25,0,0]);
deformed_tris = [ shiftmesh.v(mesh.f(:,1),:), shiftmesh.v(mesh.f(:,2),:), shiftmesh.v(mesh.f(:,3),:) ];
boundary_cons = [Bi, shiftmesh.v(Bi,:) ];


% bent mesh
deformed_tris = [ bentmesh.v(mesh.f(:,1),:), bentmesh.v(mesh.f(:,2),:), bentmesh.v(mesh.f(:,3),:) ];
boundary_cons = [Bi, bentmesh.v(Bi,:) ];


% gaussian bump
tmpmesh = mesh;
for i = 1:numel(tmpmesh.vidx)
    d2 = tmpmesh.v(i,1)^2 + tmpmesh.v(i,3)^2;
    g = 0.5 * exp(-d2/0.25);
    tmpmesh.v(i,:) = tmpmesh.v(i,:) + [0,g,0];
end
plotMesh(tmpmesh);
deformed_tris = [ tmpmesh.v(mesh.f(:,1),:), tmpmesh.v(mesh.f(:,2),:), tmpmesh.v(mesh.f(:,3),:) ];
boundary_cons = [Bi, tmpmesh.v(Bi,:) ];


deformed_mesh = poissonMesh_Yu04(mesh, boundary_cons, deformed_tris);
plotMesh(deformed_mesh);


%% curveDeform_Nealen07 test   (not implemented yet)

% make polyline
angles = ( 0:4:360 ) * pi / 180;
pl = [cos(angles'), sin(angles'), zeros(numel(angles),1)];


roi_start = 5;
roi_end = 30;
handle = 17;
displace = [0.25, 0.25, 0.0];
roi = roi_start:roi_end;
vcons = [];
vcons(1,:) = [ roi_start, pl(roi_start,:) ];
vcons(2,:) = [ roi_end,   pl(roi_end,:) ];
vcons(3,:) = [ handle,    vadd(pl(handle,:), displace) ];

pl2 = curveDeform_Nealen07(pl, roi, vcons);

newplot;
hold all;
drawpolyloop(pl, 1, 'e');
drawpolyloop(pl2, 2, 'e');
scatter3( pl(roi,1), pl(roi,2), pl(roi,3), 40 );
scatter3( vcons(:,2), vcons(:,3), vcons(:,4), 'filled' );
view([45 35])
axis equal;
hold off;
drawnow;
