%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick an input mesh to read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% portion of a sphere
mesh = readMesh('patch.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

 
mesh = readMesh('bump_ref.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

mesh = readMesh('bump_mc.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');



mesh = readMesh('patch2.obj', 'n');
plotMesh(mesh,'efb');

mesh = readMesh('bunny.obj');
plotMesh(mesh,'efb');

% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
plotMesh(mesh,'efb');

mesh = readMesh('doghead_cut.obj');
plotMesh(mesh,'efb');

mesh = readMesh('doghead_cut2.obj');
mesh = clipEars(mesh);
plotMesh(mesh,'efb');

mesh = readMesh('doghead_basehole.obj');
plotMesh(mesh,'efb');

mesh = readMesh('4bump.obj');
plotMesh(mesh,'efb');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick a weight matrix type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% one-ring weights
W = makeOneRingWeights(mesh, 'dcp');          % cotan weights
W = makeOneRingWeights(mesh, 'uniform');
W = makeOneRingWeights(mesh, 'gaussian');     % heat-kernel weights from LEM paper
W = makeOneRingWeights(mesh, 'optimal3');     % LLE weights

% epsilon-ball weights
edgelen_avg = edgeStats(mesh.v, mesh.e);
W = makeEpsBallWeights(mesh, 'uniform', 2.5*edgelen_avg);
W = makeEpsBallWeights(mesh, 'gaussian', 2.5*edgelen_avg);
W = makeEpsBallWeights(mesh, 'gaussian', 1.5*edgelen_avg);

% expmap weights (does not depend on mesh)

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 1.5*edgelen_avg;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.0*edgelen_avg;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 15;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);


% estimate boundary points from local expmaps
iboundary = zeros(size(mesh.vidx));
for i = 1:numel(mesh.vidx)
    idx = find(Vu(i,:));
    x = [0, full(Vu(i,idx))];
    y = [0, full(Vv(i,idx))];
    dists = vmag2([x',y']);
    r = 1.5*mean(dists);
   
    % voronoi(x,y);   % plot 2D voronoi diagram
    [v,c] = voronoin([x(:) y(:)]);
    v1 = v(c{1},:);
    dists = vmag2(v1);
    if max(dists) > r
        iboundary(i) = 1;
    end
end
estboundaryv = find(iboundary==1);
hold all;
plotMesh(mesh);
scatter3(mesh.v(estboundaryv,1),mesh.v(estboundaryv,2),mesh.v(estboundaryv,3));
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute LLE or LEM solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% LLE

% [RMS] increasing fudge factor in embedLLE.m will get us closer to 'flat'
%    embedding (although seems to increase distortion)
% [RMS] using boundary matrix B instead of weight-sum matrix D gives much
%    better results (need to use fudge factor as distortion gets higher)

% original LLE
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W );
subplot(1,2,1);
plotMesh(mesh, 'uefb');
subplot(1,2,2);
plotMesh(tmp,'efb');

% boundary constraints from [Mullen02]
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W, mesh.vidx(mesh.isboundaryv~=0) );
subplot(1,2,1);
plotMesh(mesh, 'uefb');
subplot(1,2,2);
plotMesh(tmp,'efb');


% boundary constraints from [Mullen02]
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W, estboundaryv );
subplot(1,2,1);
plotMesh(mesh, 'uefb');
subplot(1,2,2);
plotMesh(tmp,'efb');





% LEM

% original LEM
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W );
subplot(1,2,1);
plotMesh(mesh, 'uefb');
subplot(1,2,2);
plotMesh(tmp,'efb');

% boundary constraints from [Mullen02]
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W, mesh.vidx(mesh.isboundaryv~=0) );
subplot(1,2,1);
plotMesh(mesh, 'uefb');
subplot(1,2,2);
plotMesh(tmp,'efb');


[mesh.u,tmp.v] = embedLEM(mesh.v, W, estboundaryv );

% modified SCP with custom weight matrices

subplot(1,1,1);
mesh.u = embedSCP(mesh, 'lle', W);
plotMesh(mesh, 'uefb');

subplot(1,1,1);
mesh.u = embedSCP(mesh, 'fiedler', W);
plotMesh(mesh, 'uefb');

subplot(1,1,1);
mesh.u = embedSCP(mesh, 'generalized', W);
plotMesh(mesh, 'uefb');

subplot(1,1,1);
mesh.u = embedSCP(mesh, 'robust', W);
plotMesh(mesh, 'uefb');


% standard SCP

mesh.u = embedSCP(mesh, 'fiedler');
plotMesh(mesh, 'uefb');

mesh.u = embedSCP(mesh, 'generalized');
plotMesh(mesh, 'uefb');




writeMesh(mesh, 'temp.obj');