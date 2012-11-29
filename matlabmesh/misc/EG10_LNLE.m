


mesh = readMesh('bump_ref.obj');
mesh = clipEars(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


mesh = readMesh('bump_mc.obj');
mesh = clipEars(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


mesh = readMesh('doghead_basehole.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


mesh = readMesh('doghead_cut2.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


mesh = readMesh('bunnyhead1.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


mesh = readMesh('square.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

mesh = readMesh('circle.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

mesh = readMesh('saddle.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


mesh = readMesh('nuHemi1.obj', 'n');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');



mesh = readMesh('decimated64_T9_bp2percent.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
figure;   plotMesh(mesh,'efb');



%%  compute weight matrices


% one-ring weights
W = makeOneRingWeights(mesh, 'dcp');          % cotan weights
W = makeOneRingWeights(mesh, 'uniform');
W = makeOneRingWeights(mesh, 'gaussian');     % heat-kernel weights from LEM paper

% epsilon-ball weights
edgelen_avg = edgeStats(mesh.v, mesh.e);
W = makeEpsBallWeights(mesh.v, 'uniform', 2.5*edgelen_avg);
W = makeEpsBallWeights(mesh.v, 'gaussian', 2.5*edgelen_avg);
W = makeEpsBallWeights(mesh.v, 'gaussian', 1.5*edgelen_avg);


% DEM weights, K=10
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 15;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(mesh, options);
optimization_init_val=10;
%optimization_init_val = 75;  %?? this works for saddle...

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(mesh, options);
Wsave_dem_geoball = W;
W = Wsave_dem_geoball;  % use to restore if we computed earlier
optimization_init_val=15;



%
%  3D-optimal weights
%
% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal3';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err,Vu,Vv] = makeExpMapWeights(mesh, options);



% estimate boundary points
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,iboundary] = pointArea(P, [], 'uvVoronoi', aoptions);
estboundaryv = find(iboundary==1);
figure;
hold all;
plotMesh(mesh,'ef');
scatter3(P(estboundaryv,1),P(estboundaryv,2),P(estboundaryv,3));
hold off;
isboundaryv = zeros(N,1);
isboundaryv(estboundaryv) = 1;


loops = findBoundaries(P, estboundaryv, Vu,Vv);
figure;
hold all;
plotMesh(mesh,'ef');
for k = 1:numel(loops)
    plotPolyline(P, loops{k});
end
hold off;  drawnow;





% cotangent
W = makeOneRingWeights(mesh, 'dcp');
optimization_init_val = 1;

% normalized cotangent
W = makeOneRingWeights(mesh, 'dcp', 1);
optimization_init_val = 6;


% area-weighted cotangent
W = cotanWeights(mesh, [], 0, 1);
optimization_init_val = 1;


% uniform
W = makeOneRingWeights(mesh, 'uniform');
optimization_init_val = 1;

% normalized uniform
W = makeOneRingWeights(mesh, 'uniform', 1);
optimization_init_val = 4.4;

% cotangent-scaled uniform (gives same result as normalization)
W = makeOneRingWeights(mesh, 'uniform');
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    if isboundaryv(k)
        valence = 6;
    end
    angle = (pi/2 - pi/valence);
    W(k,idx) = 2 * cot(angle);
end
optimization_init_val = 1;

% mean-value 
W = makeOneRingWeights(mesh, 'meanvalue');
optimization_init_val = 1;

% mean-value (normalized)
W = makeOneRingWeights(mesh, 'meanvalue', 1);
optimization_init_val = 1;

% authalic
W = makeOneRingWeights(mesh, 'dap');
optimization_init_val = 1;

% eps-ball uniform
edgelen_avg = edgeStats(mesh.v, mesh.e);
W = makeEpsBallWeights(mesh.v, 'uniform', 1.5*edgelen_avg);

% heat-kernel one-ring weights 
W = makeOneRingWeights(mesh, 'gaussian');     

% LLE k-nbr weights
W = makeKNbrWeights(mesh.v, 'optimal3', 8);


% check that weights sum to 1, and check reconstruction error
D = -W;
D(1:N+1:N*N) = -sum(D,2);
sum(sum(D,2))
sum(D*mesh.v)


% compute valence from weight matrix
Wone = W; k=find(Wone(:));  Wone(k) = Wone(k) ./ Wone(k);
valence = full(sum( Wone, 2 ));

% scale weights by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    W(k,idx) = valence * W(k,idx);
end
% scale interior weights
for k = 1:N
    idx = find(W(k,:));
    if ( ~ isboundaryv(k) )
        W(k,idx) = 8 * W(k,idx);
    end
end
%scale interior by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    if ( ~ isboundaryv(k) )
        W(k,idx) = valence * W(k,idx);
    end
end


% do not use specific constraint points
constrain_vtx1 = -1;
constrain_vtx2 = -1;


% compute solution just to figure out which points to constrain
[mesh.u,EC,cons1,cons2] = embedDNCP(mesh, W);
figure; newplot; hold on; 
plotMesh(mesh, 'uefb');
consv = mesh.u([cons1,cons2],:);
scatter(consv(:,1), consv(:,2), 200, 'ms', 'filled');
hold off;  drawnow;
constrain_vtx1 = cons1;
constrain_vtx2 = cons2;


% plot cotan DNCP as we vary which pair of boundary points are constrained
figure;
Wcot = cotanWeights(mesh);
for k = 2:numel(mesh.loops)
    [mesh.u,Ed,cons1,cons2] = embedDNCP(mesh, Wcot, mesh.loops(1), mesh.loops(k));
    newplot; hold on; 
    plotMesh(mesh, 'uefb');
    consv = mesh.u([cons1,cons2],:);
    scatter(consv(:,1), consv(:,2), 200, 'ms', 'filled');
    hold off;
    drawnow;
    pause(0.1);
end




% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
figure; plotMesh(mesh, 'UefbO');




% find optimal W scale value k  using SCP
init_val = optimization_init_val;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, x, W ), init_val/2, 2*init_val );
mesh.u = embedSCP(mesh, 'generalized', k*W );
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
if exist('myjet','var')  colormap(myjet); caxis([cmin,cmax]);  end
title(['Optimized SCP (QC:', num2str(sum(faceQC)),')']);



% find optimal W scale value k  using DNCP
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W), 0.001, 10 );
mesh.u = embedDNCP(mesh, k*W );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
figure; plotMesh(mesh, 'Uefb');




X = [];
Cvals = [];    Dvals = [];   Avals = [];
kCvals = [];   kDvals = [];  kAvals = [];
Evals = [];
turnVals = [];
%for k = 0.9:0.01:1.5       % uniform weights
%for k = 4:0.01:8
%for k = 0.0777:0.00001:0.78
%for k = 1:0.1:10
%for k = 3:0.01:5           % normalized uniform weights (broad)
%for k = 4.0:0.001:4.2           % normalized cotan weights (fine)
%for k = 0.9:0.001:1.2   % show turning number jumps for cotan weights
%for k = 0.86:0.0001:1.1
%for k = 0.975:0.00001:1.025
for k = 0.998:0.00001:1.00
%for k = 0.997:0.00001:1.001   % show turning number jumps for cotan weights
%for k = 0.9985:0.00001:1.0     % show problem for saddle

    mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2);
    E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
    
    %mesh.u = embedSCP(mesh, 'generalized', k*W);
    %E = FnOptLaplacianScale( 'scp', mesh, k, W, constrain_vtx1, constrain_vtx2 );

    turnVals = [turnVals, turningNumber(mesh.u(mesh.loops,:))];
    [c,d,a] = eConformal(mesh, W);
    [c,d,a,t] = eConformal(mesh, k*W);

    X = [X,k];
    Cvals = [Cvals,c];   Dvals = [Dvals,d];  Avals = [Avals,a];
    kCvals = [kCvals,c];   kDvals = [kDvals,d];  kAvals = [kAvals,a];
    Evals = [Evals,E];
    
    plotMesh(mesh, 'Uefb');    
%    xlim([-2,2]);
    drawnow; %pause;
end

Ascale = kAvals.^-1;

figure;
hold all;
Cscale = mean(kCvals.^2);
plot(X,Cscale*turnVals);
plot(X,kCvals.^2);
plot(X,kCvals);
plot(X,abs(kCvals));
hold off;
ylim([-0.001,0.001])
drawnow;


figure;
hold all;
plot(X,Dvals);
plot(X,Avals);
plot(X,Cvals);
Cscale = mean(plotvals)
plot(X,Cscale*turnVals);
hold off;
drawnow;


k=0.9995;
k = 44;
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2);
plotMesh(mesh, 'Uefb');    

    
% [SADDLE]
mesh.u = embedDNCP(mesh, 22*W, 1);
figure; plotMesh(mesh, 'uefb');
    

mesh.u = embedSCP(mesh, 'generalized', 1.5*W);
figure; plotMesh(mesh, 'Uefb');









%%  SCP embedding for meshes

mesh.u = embedSCP(mesh, 'fiedler');
figure;     plotMesh(mesh, 'Uefb');
title('SCP - Fiedler Vector');

mesh.u = embedSCP(mesh, 'generalized');
plotMesh(mesh, 'UefbO');
%figure;     plotMesh(mesh, 'UefbO');
title('SCP - Generalized Eigen');

mesh.u = embedDNCP(mesh);
figure;     plotMesh(mesh, 'Uefb');


% area-weighted SCP
W = cotanWeights(mesh, [], 0, 1);
mesh.u = embedSCP(mesh, 'generalized',W, 1);
figure;     plotMesh(mesh, 'Uefb');
title('SCP - Generalized Eigen');

%dogface
mesh.u = embedSCP(mesh, 'generalized',13.906907*W);
figure;     plotMesh(mesh, 'UefbO');
title('SCP - Generalized Eigen');

tmp = mesh;
tmp.v = [ mesh.u(:,1), mesh.u(:,2), zeros(N,1) ];
A = faceArea(tmp);

meshSCP = mesh;
meshSCP.u = embedSCP(meshSCP, 'generalized');
figure; plotMesh( meshSCP, 'Uefb');
meshDNCP = mesh;
%Wcot = cotanWeights(mesh);
Wcot = cotanWeights(mesh) * 0.999529;
meshDNCP.u = embedDNCP(mesh, Wcot);
figure; plotMesh( meshDNCP, 'Uefb');
[EC_SCP,ED_SCP,A_SCP] = eConformal(meshSCP,Wcot)
[EC_DNCP,ED_DNCP,A_DNCP] = eConformal(meshDNCP,Wcot)
fprintf('[Spectral] C: %5.5f   D: %5.5f  A: %5.5f\n', EC_SCP, ED_SCP, A_SCP);
fprintf('[Natural ] C: %5.5f   D: %5.5f  A: %5.5f\n', EC_DNCP, ED_DNCP, A_DNCP);


% hack - SCP w/ custom weights. By scaling weights
%  we can get quite nice parameterizations, even though 
%  the resulting ((D-W)-A) matrix is not positive definite (weird...)
Wscale = 1/mean(nonzeros(W(:)));
figure;
mesh.u = embedSCP(mesh, 'generalized', 16*W);
plotMesh(mesh, 'uefb');




%% LLE

% [RMS] increasing fudge factor in embedLLE.m will get us closer to 'flat'
%    embedding (although seems to increase distortion)
% [RMS] using boundary matrix B instead of weight-sum matrix D gives much
%    better results (need to use fudge factor as distortion gets higher)

% original LLE
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W );
figure;
subplot(1,2,1);
plotMesh(mesh, 'Uefb');
title('LLE - 2D');
subplot(1,2,2);
plotMesh(tmp,'efb');
title('LLE - 3D');

% original LLE (no 3D plot)
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W );
figure;
plotMesh(mesh, 'Uefb');
title('LLE - 2D');

% boundary constraints from [Mullen02]  (use real boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W, mesh.vidx(mesh.isboundaryv~=0) );
figure;
plotMesh(mesh, 'Uefb');
title('LLE+RealBoundary - 2D');


% boundary constraints from [Mullen02]   (use estimated boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W, estboundaryv );
figure;
plotMesh(mesh, 'Uefb');
title('LLE+EstBoundary - 2D');


% mesh boundary points + loop
mesh.u = embedLLE(mesh.v, W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
figure;
plotMesh(mesh, 'Uefb');
title('LLE+EstBoundary - 2D');



% C-LLE - optimize w/ mesh boundary 
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'lle', mesh, x, W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh ), init_val/2, 2*init_val );
mesh.u = embedLLE(mesh.v, k*W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)),')']);



% C-LLE - optimize w/ estimated boundary points + estimated loop
TRIS = delaunay(mesh.u(:,1), mesh.u(:,2));
delmesh = makeMesh(P, TRIS);

% compute initial param
mesh.u = embedLLE(mesh.v, W, delmesh.vidx(delmesh.isboundaryv~=0), delmesh.loops );
figure;
plotMesh(mesh, 'Uefb');
title('C-LLE+EstBoundary - 2D');

% optimize
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'lle', mesh, W, x, estboundaryv, delmesh.loops ), init_val/2, 2*init_val );
mesh.u = embedLLE(mesh.v, k*W, estboundaryv, delmesh.loops );
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)),')']);





Cvals = [];    Dvals = [];   Avals = [];
kCvals = [];   kDvals = [];  kAvals = [];
X = [];  Evals = [];   turnVals = [];
for k = 5.4:0.1:6
    mesh.u = embedLLE(mesh.v, k*W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
    E = FnOptLaplacianScale( 'lle', mesh, k, W, mesh.vidx(mesh.isboundaryv~=0) );

    turnVals = [turnVals, turningNumber(mesh.u(mesh.loops,:))];
    [c,d,a] = eConformal(mesh, W);
    [c,d,a,t] = eConformal(mesh, k*W);
    X = [X,k];   Evals = [Evals,E];
    Cvals = [Cvals,c];   Dvals = [Dvals,d];  Avals = [Avals,a];
    kCvals = [kCvals,c];   kDvals = [kDvals,d];  kAvals = [kAvals,a];
    
    plotMesh(mesh, 'Uefb');      drawnow;
end


plotvals = kCvals.^2;
figure;
hold all;
Cscale = mean(plotvals)
plot(X,Cscale*turnVals);
plot(X,plotvals);
hold off;
drawnow;






init_val = 1/mean(nonzeros(W(:))) * 1.7;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'lle', mesh, W, x ), init_val/2, 2*init_val );
mesh.u = embedLLE(mesh.v, k*W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
fprintf('k: %f   C: %f   turning: %f\n', k, EConformal, Turning);
figure; plotMesh(mesh, 'Uefb');
    




%% LEM

% original LEM
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W );
figure;
plotMesh(mesh, 'Uefb');
title('LEM - 2D');



% boundary constraints from [Mullen02]  (use real boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W, mesh.vidx(mesh.isboundaryv~=0) );
figure;
plotMesh(mesh, 'Uefb');
title('LEM+RealBoundary - 2D');


% boundary constraints from [Mullen02]  (use estimated boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W, estboundaryv );
figure;
plotMesh(mesh, 'Uefb');
title('LEM+RealBoundary - 2D');



%% LTSA

mesh.u = ltsa(P', 2, 8)';
figure;
plotMesh(mesh, 'UefbO');

mesh.u = HessianLLE(P', 8, 2)';
figure;
plotMesh(mesh, 'UefbO');

mesh.u = HessianLLE(P', 8, 2, W')';
figure;
plotMesh(mesh, 'UefbO');

mesh.u = LEigenmaps(P, 8, 2);
figure;
plotMesh(mesh, 'UefbO');


edgelen_avg = edgeStats(mesh.v, mesh.e);
mesh.u = compute_mapping(P, 'GPLVM', 2, edgelen_avg*5);
figure;
plotMesh(mesh, 'UefbO');

edgelen_avg = edgeStats(mesh.v, mesh.e);
mesh.u = compute_mapping(P, 'tSNE', 2);
figure;
plotMesh(mesh, 'UefbO');





%% remesh based on parameterization

TRIS = delaunay(mesh.u(:,1), mesh.u(:,2));
tmp = makeMesh(P, TRIS);
tmp.u = mesh.u;
figure;
plotMesh(tmp, 'efb');
title('3D remesh');


figure;
hold all;
plotMesh(tmp, 'uef');
scatter(mesh.u(:,1),mesh.u(:,2));
title('UV Remesh');
hold off;


tmp.isboundaryv


voronoi(mesh.u(:,1),mesh.u(:,2));
%[V,C]  = voronoin(mesh.u);




%% make non-uniform hemisphere mesh


    
uv1 = genpoints(1000,2,'stratified_unit_square');    
uv1 = vadd(uv1, [-0.5,-0.5])*2;
vdists = vmag2(uv1);
uv1 = uv1(vdists<1,:);
figure; scatter(uv1(:,1),uv1(:,2));  axis equal;

uv2 = genpoints(4000,2,'stratified_unit_square');    
uv2 = vadd(uv2, [-0.5,-0.5])*2;
vdists = vmag2(uv2);
uv2 = uv2(vdists<1 & uv2(:,1) > 0 ,:);
figure; scatter(uv2(:,1),uv2(:,2));  axis equal;

uv = [uv1;uv2];

points = [uv(:,1), uv(:,2), sqrt(ones(size(uv,1),1) - vmag2(uv))];
figure; scatter3(points(:,1),points(:,2),points(:,3));  axis equal;

TRIS = delaunay(uv(:,1), uv(:,2));
spmesh = makeMesh(points, TRIS);
plotMesh(spmesh);


writeMesh(spmesh,'nuHemi1.obj');





%% highly non-uniform point distributions

ps1 = readPoints('tmp_sphere1.obj');
ps2 = readPoints('tmp_sphere2.obj');

idx1 = find(ps1.v(:,3) > 0 & ps1.v(:,1) > 0);
idx2 = find(ps2.v(:,3) > 0 & ps2.v(:,1) < 0);

pts = normalize([ps1.v(idx1,:); ps2.v(idx2,:)]);
ps = makePointSet(pts, pts);
plotPoints(ps, 'v');

edgelen_avg = edgeStats(ps.v, ps.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3*edgelen_avg;
W = makeExpMapWeights(ps, options);

ps.u = embedLLE(ps.v, W);
plotPoints(ps,'uv');

% generate topology
TRIS = delaunay(uv(:,1),uv(:,2));
mesh = makeMesh(ps.v, TRIS, ps.n);
plotMesh(mesh, 'efb');


mesh.u = embedSCP(mesh, 'generalized', [], 0);
figure; plotMesh(mesh,'UefbO'); title('Cotan-SCP');

WcotA = cotanWeights(mesh, [], 0, 1);
mesh.u = embedSCP(mesh, 'generalized', WcotA, 1);
figure; plotMesh(mesh,'UefbO'); title('Area-Weighted Cotan-SCP');


options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 30;
[W,err,Vu,Vv] = makeExpMapWeights(ps, options);   Wsave=W;  W = Wsave;
init_val = 7;

options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3*edgelen_avg;
[W,err,Vu,Vv] = makeExpMapWeights(ps, options);    Wsave=W;  W = Wsave;
init_val = 8;

% scale rows by point area
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
A = pointArea(ps.v, [], 'uvVoronoi', aoptions);
A = vertexArea(mesh, [], 'onering');
A = A / mean(A);
for k = 1:size(W,1);
    if ( mesh.isboundaryv(k) )
        continue;
    end
    idx = find(W(k,:));
%    scale = numel(idx);
    scale = A(k);
    W(k,idx) = scale * W(k,idx);
end
init_val = 10;

% scale entries by point area
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
A = pointArea(ps.v, [], 'uvDelRing', aoptions);
A = A / mean(A);
for k = 1:size(W,1);
    idx = find(W(k,:));
    for j = 1:numel(idx)
        W(k,idx(j)) = A(idx(j)) * W(k,idx(j));
    end
end
init_val = 10;


% find optimal W scale value k  using SCP
ps.loops{1} = mesh.loops{1};
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'pscp', ps, x, W ), init_val/4, 4*init_val );
mesh.u = embedPSCP(ps, k*W );
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
if exist('myjet','var')  colormap(myjet); caxis([cmin,cmax]);  end
title(['Optimized PSCP (QC:', num2str(sum(faceQC)),')']);

plotMesh(mesh,'efbO',A);

ps.u = embedLLE(ps.v, W);
plotPoints(ps,'uv');

