
pointset = readPoints('dogface.pts');


pointset = readPoints('vertebra_top.pts');


pointset = readPoints('vertebra_top_N.pts');
N = size(pointset.v,1);


pointset = readPoints('face_5k.pobj');
N = size(pointset.v,1);



lightdir = normalize([0.5,0.5,0.5]);
NdotL = abs(vdot(pointset.n, lightdir));
ptColors = [NdotL,NdotL,NdotL];
scatter3(pointset.v(:,1),  pointset.v(:,2), pointset.v(:,3), 25, ptColors, 'Marker','.');
axis equal;



% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(pointset.v, pointset.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);


% DEM weights, K=15
edgelen_avg = edgeStats(pointset.v, pointset.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 30;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);


fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'k', 30);

fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 4*8.3);


% estimate boundary points
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);


% figure;
% hold on;
% for i = 1:numel(estboundaryv)
%     vi = estboundaryv(i); 
%     idx = find(Vu(vi,:));
%     x = [0, full(Vu(vi,idx))];
%     y = [0, full(Vv(vi,idx))];
%     uv = [x',y'];
%     TRI = delaunay(x,y);
%     %xyz = [uv, zeros(size(x))'];
%     xyz = [pointset.v(vi,:); pointset.v(idx,:)];
%     tmpmesh = makeMesh(xyz,TRI);
%     plotMesh(tmpmesh, 'ef');
%     scatter3(pointset.v(vi,1),pointset.v(vi,2),pointset.v(vi,3));
% end
% hold off;
% drawnow;

angle_gaps = zeros(N,1);
for i = 1:N
    idx = find(Vu(i,:));
    x = [0, full(Vu(i,idx))];
    y = [0, full(Vv(i,idx))];
    uv = [x',y'];
    angles = atan2(y,x) + pi;
    [angles,sidx] = sort(angles);
    maxgap = 0;
    for k = 2:numel(angles)
        gap = angles(k) - angles(k-1);
        maxgap = max(maxgap,gap);
    end
    maxgap = max(maxgap, 2*pi+angles(1) - angles(end));
    angle_gaps(i) = maxgap;
end

scatter3(pointset.v(:,1),  pointset.v(:,2), pointset.v(:,3), 25, angle_gaps,'Marker','.');
colorbar;
axis equal;

estboundaryv = find(angle_gaps > pi*0.65);

figure;
hold all;
lightdir = normalize([0.5,0.5,0.5]);
NdotL = abs(vdot(pointset.n, lightdir));
ptColors = [NdotL,NdotL,NdotL];
%scatter3(pointset.v(:,1),  pointset.v(:,2), pointset.v(:,3), 25, ptColors, 'Marker','.');
scatter3(pointset.v(estboundaryv,1),pointset.v(estboundaryv,2),pointset.v(estboundaryv,3));
axis equal;
hold off;



loops = findBoundaries(pointset, estboundaryv, '3D', Vu,Vv);
pointset.loops = {loops{1}};
figure;
hold all;
%plotMesh(mesh,'ef');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
axis equal;
hold off;  drawnow;




%% face 40k


pointset = readPoints('face_40k.pobj');
N = size(pointset.v,1);

% estimate boundary points
edgelen_avg = edgeStats(pointset.v, pointset.e);
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 2.1*edgelen_avg);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);

loops = findBoundaries(pointset, estboundaryv, '3D');
pointset.loops = {loops{1}};
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
axis equal;  hold off;  drawnow;


% compute weights

% DEM weights, G=2.1*edgelenavg
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);



% original LEM
pointset.u = embedLEM(pointset.v, W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n);
figure; plotMesh(delmesh,'efbl');  title('LEM mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM uv');


% LEM with boundary
pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');


% optimized LEM
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 10,30 );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');





%% face 10k

pointset = readPoints('face_10k.pobj');
N = size(pointset.v,1);

% estimate boundary points
edgelen_avg = edgeStats(pointset.v, pointset.e);
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 2.5*edgelen_avg);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);

loops = findBoundaries(pointset, estboundaryv, '3D');
pointset.loops = {loops{1}};
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
axis equal;  hold off;  drawnow;


% compute weights

% DEM weights, G=2.1*edgelenavg
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);



% original LEM
pointset.u = embedLEM(pointset.v, W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM uv');
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;

% writeMesh(delmesh, 'face10k_LEM.obj');


% LEM with boundary
pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;

% writeMesh(delmesh, 'face10k_LEMB.obj');


% optimized LEM
options = struct('TolX',sqrt(eps));
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 10,30, options );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;


% writeMesh(delmesh, 'face10k_PSCP.obj');


mesh = readMesh('face_10k_bp2p5.obj');
mesh.u = embedSCP(mesh);
figure; plotMesh(mesh,'ufb', mesh.n);   title('SCP uv');



%% face 5K


pointset = readPoints('face_5k.pobj');
N = size(pointset.v,1);


% estimate boundary points
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 4*8.3);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);

loops = findBoundaries(pointset, estboundaryv, '3D');
pointset.loops = {loops{1}};
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
axis equal;  hold off;  drawnow;


% compute weights

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(pointset.v, pointset.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.1*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);



% boundary LEM


pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;

% writeMesh(delmesh, 'face5k_LEMB.obj');



[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 20, 30 );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;

% writeMesh(delmesh, 'face5k_PSCP.obj');


TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS);
delmesh.n = estimateNormal(delmesh);
figure; plotMesh(delmesh,'efbl');








%% gnome 13k

pointset = readPoints('gnome_13k.pobj');
N = size(pointset.v,1);

% estimate boundary points
edgelen_avg = edgeStats(pointset.v, pointset.e);
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 2.5*edgelen_avg);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);

loops = findBoundaries(pointset, estboundaryv, '3D');
pointset.loops = {loops{1}};
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
axis equal;  hold off;  drawnow;

load('myjet.mat');  cmin = 1;   cmax = 1.5;

% compute weights

% DEM weights, G=2.1*edgelenavg
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);



% original LEM
pointset.u = embedLEM(pointset.v, W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('LEM mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmeshuv,'f', delmesh.n);   title('LEM uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
hold off; drawnow;
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;

% writeMesh(delmesh, 'gnome_LEM.obj');

% LEM with boundary
pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('LEM+estboundary mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmeshuv,'f', delmesh.n);   title('LEM+estboundary uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
hold off; drawnow;
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;


% writeMesh(delmesh, 'gnome_LEMB.obj');


% optimized LEM
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 10,30 );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('LEM+estboundary mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmeshuv,'f', delmesh.n);   title('LEM+estboundary uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
hold off; drawnow;
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;


% writeMesh(delmesh, 'gnome_PSCP.obj');


% clip tris on flattened mesh
keepidx = [];
xv = delmesh.u(pointset.loops{1},1);
yv = delmesh.u(pointset.loops{1},2);
xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
for ti = 1:numel(delmesh.fidx)
    centroid = sum(delmesh.u(delmesh.f(ti,:),:)) / 3;
    if inpolygon(centroid(1),centroid(2), xv, yv)
        keepidx = [keepidx;ti];
    end
end
faces = delmesh.f(keepidx,:);
%faces = [faces(:,2), faces(:,1), faces(:,3)];
delmesh2 = makeMesh(delmesh.v, faces, delmesh.n, delmesh.u);
plotMesh(delmesh2);

delmesh2.u = embedSCP(delmesh2, 'generalized', [], 0);
plotMesh(delmesh2,'uefb');


err = faceDistortion(delmesh2, 'qc');
hold off; drawnow;
figure; hold on;
plotMesh(delmesh2, 'uef', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmesh2.u, pointset.loops{k});
end
hold off; drawnow;



%% igea 11k

pointset = readPoints('igea_11k.pobj');
N = size(pointset.v,1);

radius = 150;
ballpt = [0,0,-290];
dists = vmag2(vadd(pointset.v,-ballpt));
usepts = find(dists>radius*radius);
usev = pointset.v(usepts,:);
usen = pointset.n(usepts,:);
pointset = makePointSet(usev,usen);
N = size(pointset.v,1);

% estimate boundary points
edgelen_avg = edgeStats(pointset.v, pointset.e);
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 3.5*edgelen_avg);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);
estboundaryv = estboundaryv(estboundaryv~=2441);   % hack: fix problematic vert...


loops = findBoundaries(pointset, estboundaryv, '3D', [],[], 7);
pointset.loops = {loops{1}};
%pointset.loops = loops;
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
scatter3(pointset.v(estboundaryv,1),pointset.v(estboundaryv,2),pointset.v(estboundaryv,3));
axis equal;  hold off;  drawnow;


% compute weights

% DEM weights, G=2.1*edgelenavg
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.5*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);

load('myjet.mat');  cmin = 1;   cmax = 1.5;

% original LEM
pointset.u = embedLEM(pointset.v, W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('LEM+estboundary mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmesh, 'uf', delmesh.n);   title('LEM+estboundary uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;




% LEM with boundary
pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('LEM+estboundary mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmesh, 'uf', delmesh.n);   title('LEM+estboundary uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;



% writeMesh(delmesh, 'igea_LEMB.obj');



% optimized LEM
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 2,5 );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('Wtan-SCP mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmesh, 'uf', abs(delmesh.n));   title('Wtan-SCP uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uef', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;



% writeMesh(delmesh, 'igea_PSCP.obj');










%% igea 33k

pointset = readPoints('igea_33k.pobj');
N = size(pointset.v,1);

radius = 150;
ballpt = [0,0,-290];
dists = vmag2(vadd(pointset.v,-ballpt));
usepts = find(dists>radius*radius);
usev = pointset.v(usepts,:);
usen = pointset.n(usepts,:);
pointset = makePointSet(usev,usen);
N = size(pointset.v,1);

% estimate boundary points
edgelen_avg = edgeStats(pointset.v, pointset.e);
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 3.5*edgelen_avg);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);
estboundaryv = estboundaryv(estboundaryv~=12758);   % hack: fix problematic vert...
estboundaryv = estboundaryv(estboundaryv~=12650);   % hack: fix problematic vert...
estboundaryv = estboundaryv(estboundaryv~=23173);   % hack: fix problematic vert...


loops = findBoundaries(pointset, estboundaryv, '3D', [],[], 6);
pointset.loops = {loops{1}};
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
text( pointset.v(estboundaryv,1),pointset.v(estboundaryv,2),pointset.v(estboundaryv,3), int2str(estboundaryv) );
scatter3(pointset.v(estboundaryv,1),pointset.v(estboundaryv,2),pointset.v(estboundaryv,3));
axis equal;  hold off;  drawnow;


% compute weights

% DEM weights, G=2.1*edgelenavg
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);

load('myjet.mat');  cmin = 1;   cmax = 1.5;




% LEM with boundary
pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('LEM+estboundary mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmesh, 'uf', abs(delmesh.n));   title('LEM+estboundary uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;


% writeMesh(delmesh, 'igea33k_LEMB.obj');





% optimized LEM
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 10,30 );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure;  hold on;
plotMesh(delmesh,'ef');  title('Wtan-SCP mesh');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off; drawnow;
figure; hold on;
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
plotMesh(delmesh, 'uf', abs(delmesh.n));   title('Wtan-SCP uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uef', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;



% writeMesh(delmesh, 'igea33k_PSCP.obj');


%% bump boundary example

% colormap info
cmin = 1;   cmax = 2;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');


mesh = readMesh('bump_ref.obj');
mesh = clipEars(mesh);
meshboundaryv = mesh.loops{1};
pointset = makePointSet(mesh.v,mesh.n);
N = size(pointset.v,1);
plotMesh(mesh);

edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);


% estimate boundary points
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);
%figure;
%hold all;
%plotMesh(mesh,'ef');
%scatter3(pointset.v(estboundaryv,1),pointset.v(estboundaryv,2),pointset.v(estboundaryv,3));
%hold off;


pointset.loops = findBoundaries(pointset.v, estboundaryv, Vu,Vv);
figure;
hold all;
plotMesh(mesh,'ef');
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
hold off;  drawnow;



mesh.u = embedLEM(mesh.v, W, meshboundaryv );
figure;
faceQC = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title('LEM / mesh boundary');


mesh.u = embedLEM(mesh.v, W, estboundaryv );
figure;
hold on;
faceQC = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'uefO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline([mesh.u,zeros(N,1)], pointset.loops{k});
end
title('LEM / estimated boundary');
hold off;


tmp = pointset;
tmp.loops{1} = pointset.loops{1}(1:2:end);
mesh.u = embedLEM(mesh.v, W, tmp.loops{1} );
figure;
hold on;
faceQC = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'uefO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(tmp.loops)
    plotPolyline([mesh.u,zeros(N,1)], tmp.loops{k});
end
title('LEM / reduced estimated boundary');
hold off;



% use angle_gap technique

angle_gaps = zeros(N,1);
for i = 1:N
    idx = find(Vu(i,:));
    x = [0, full(Vu(i,idx))];
    y = [0, full(Vv(i,idx))];
    uv = [x',y'];
    angles = atan2(y,x) + pi;
    [angles,sidx] = sort(angles);
    maxgap = 0;
    for k = 2:numel(angles)
        gap = angles(k) - angles(k-1);
        maxgap = max(maxgap,gap);
    end
    maxgap = max(maxgap, 2*pi+angles(1) - angles(end));
    angle_gaps(i) = maxgap;
end
estboundaryv = find(angle_gaps > pi*0.25);

figure;
hold all;
plotMesh(mesh,'ef');
scatter3(pointset.v(estboundaryv,1),pointset.v(estboundaryv,2),pointset.v(estboundaryv,3));
hold off;


mesh.u = embedLEM(mesh.v, W, estboundaryv );
figure;
hold on;
faceQC = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'uefO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline([mesh.u,zeros(N,1)], pointset.loops{k});
end
title('LEM / estimated boundary');
hold off;


estboundaryv2 = estboundaryv(1:2:end);
mesh.u = embedLEM(mesh.v, W, estboundaryv2 );
figure;
hold on;
faceQC = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'uefO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(tmp.loops)
    plotPolyline([mesh.u,zeros(N,1)], tmp.loops{k});
end
title('LEM / reduced estimated boundary');
hold off;


% find optimal W scale value k  using SCP

tmp = pointset;
tmp.loops{1} = mesh.loops{1};
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', tmp, x, W ), 15, 25 );
mesh.u = embedPSCP(tmp, k*W );
figure;
plotMesh(mesh,'uefbO');
title('PSCP / mesh boundary');


tmp = pointset;
tmp.loops{1} = pointset.loops{1};
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', tmp, x, W ), 5, 25 );
mesh.u = embedPSCP(tmp, k*W );
figure;
hold on;
plotMesh(mesh,'uefO');
for k = 1:numel(pointset.loops)
    plotPolyline([mesh.u,zeros(N,1)], pointset.loops{k});
end
title('PSCP / estimated boundary');
hold off;






mesh.u = embedSCP(mesh,'generalized');
figure;
hold on;
faceQC = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title('mesh SCP');






%% PLOT FUNCTIONS

cmin = 1;   cmax = 1.5;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');

figure;
err = faceDistortion(mesh, 'qc');
plotMesh(mesh, 'UfbO', err );
colormap(myjet); caxis([cmin,cmax]);







%% face 10k - hole

mesh = readMesh('face_10k_holes.pobj');
plotMesh(mesh,'efb');
mesh.u = embedSCP(mesh);
plotMesh(mesh,'uefb');

% writeMesh(mesh, 'face10k_holes_SCP.obj');

pointset = makePointSet(mesh.v, mesh.n);
N = size(pointset.v,1);

% estimate boundary points
edgelen_avg = edgeStats(pointset.v, pointset.e);
fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);
[Vu,Vv] = fastExpMaps(fake_mesh, 'g', 2.5*edgelen_avg);

aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(pointset.v, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);

loops = findBoundaries(pointset, estboundaryv, '3D');
pointset.loops = loops;
%pointset.loops = {loops{1}};
figure;  hold all;
for k = 1:numel(pointset.loops)
    plotPolyline(pointset.v, pointset.loops{k});
end
axis equal;  hold off;  drawnow;


%pointset.loops = mesh.loops;

% compute weights

% DEM weights, G=2.1*edgelenavg
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(pointset, options);


% LEM with boundary
pointset.u = embedLEM(pointset.v, W, pointset.loops{1} );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; plotMesh(delmeshuv,'fb', delmesh.n);   title('LEM+estboundary uv');
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uf', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;

% writeMesh(delmesh, 'face10k_LEMB.obj');


% optimized LEM
options = struct('TolX',sqrt(eps));
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', pointset, x, W ), 10,30, options );
pointset.u = embedPSCP(pointset, k*W );

TRIS = delaunay(pointset.u(:,1), pointset.u(:,2));
delmesh = makeMesh(pointset.v, TRIS, pointset.n, pointset.u);
figure; plotMesh(delmesh,'efbl');  title('LEM+estboundary mesh');
delmeshuv = delmesh;
delmeshuv.v = [pointset.u,zeros(N,1)];
figure; hold on;
plotMesh(delmeshuv,'fb', abs(delmesh.n));   title('LEM+estboundary uv');
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;
err = faceDistortion(delmesh, 'qc');
figure; hold on;
plotMesh(delmesh, 'uef', err );
colormap(myjet); caxis([cmin,cmax]);
for k = 1:numel(pointset.loops)
    plotPolyline(delmeshuv.v, pointset.loops{k});
end
hold off; drawnow;


plotMesh(delmesh, 'uefi', err );
colormap(myjet); caxis([cmin,cmax]);

W2 = W;
%W2(W<0) = 1;
W2(W>0) = 0;
hold on;
plotMesh(delmesh, 'uef', abs(sum(W2,2)));
usev = delmesh.u;
for k = 1:numel(pointset.loops)
    loop = pointset.loops{k};
    plotPolyline(usev, loop);
    scatter(usev(loop,1), usev(loop,2), 50);
end
hold off;


% writeMesh(delmesh, 'face10k_holes_PSCP.obj');

v=7985;
idx = find(W(v,:));
nbrs = [idx,v];
uv = [ full(Vu(v,nbrs)); full(Vv(v,nbrs)) ]';
TRIS=delaunay(uv(:,1),uv(:,2));
triplot(TRIS,uv(:,1),uv(:,2));

