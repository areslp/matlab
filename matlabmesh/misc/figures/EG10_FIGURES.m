
%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of LEM and LLE, with and without boundary weight
%%%%%%%%%%%%%%%%%%%%%%%%

% colormap info
cmin = 1;   cmax = 1.5;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');
render_aa = 0;

% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

Wcot = makeOneRingWeights(mesh, 'dcp');

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh, options);
Wsave_dem_geoball = W;
W = Wsave_dem_geoball;  % use to restore if we computed earlier

% DEM weights, K=15
%options = [];
%options.weightmode = 'optimal2';
%options.nbrtype = 'k';
%options.nbrsize = 15;
%options.regtol = 10e-3;
%[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
%Wsave_dem_geoball = W;
%W = Wsave_dem_geoball;  % use to restore if we computed earlier

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
estboundaryv = find(angle_gaps > pi*0.65);
figure;
hold all;
plotMesh(mesh,'ef');
scatter3(P(estboundaryv,1),P(estboundaryv,2),P(estboundaryv,3));
hold off;
isboundaryv = zeros(N,1);
isboundaryv(estboundaryv) = 1;


% original LLE (no 3D plot)
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('LLE - 2D');
if render_aa myaa([4,2]); end

% boundary constraints from [Mullen02]   (use estimated boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W, estboundaryv );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('LLE+EstBoundary - 2D');
if render_aa myaa([4,2]); end

% original LEM
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('LEM - 2D');
if render_aa myaa([4,2]); end

% boundary constraints from [Mullen02]  (use estimated boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W, estboundaryv );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('LEM+RealBoundary - 2D');
if render_aa myaa([4,2]); end

%scale interior by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    if ( ~ isboundaryv(k) )
        W(k,idx) = valence * W(k,idx);
    end
end


% boundary constraints from [Mullen02]   (use estimated boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, W, estboundaryv );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('Scaled LLE+EstBoundary - 2D');
if render_aa myaa([4,2]); end

% boundary constraints from [Mullen02]  (use estimated boundary)
tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, W, estboundaryv );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('Scaled LEM+RealBoundary - 2D');
if render_aa myaa([4,2]); end



% cotan weights
mesh.u = embedSCP(mesh, 'generalized');
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('SCP - Generalized Eigen');
if render_aa myaa([4,2]); end

tmp = mesh;
[mesh.u,tmp.v] = embedLLE(mesh.v, Wcot, mesh.loops{1} );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('Cotan LLE+EstBoundary');
if render_aa myaa([4,2]); end

tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, Wcot, mesh.loops{1} );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('Cotan LEM+EstBoundary');
if render_aa myaa([4,2]); end

% uniform weights

Wuni = makeOneRingWeights(mesh, 'uniform');

[mesh.u,tmp.v] = embedLLE(mesh.v, Wuni, mesh.loops );
figure;
faceColor = faceColors(mesh, faceDistortion(mesh, 'qc'), cmin, cmax, clo, chi );
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('Uniform LLE+EstBoundary');
if render_aa myaa([4,2]); end

tmp = mesh;
[mesh.u,tmp.v] = embedLEM(mesh.v, Wuni, mesh.loops );
figure;
faceColor = faceColors(mesh, faceDistortion(mesh, 'qc'), cmin, cmax, clo, chi );
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('Uniform LEM+EstBoundary');
if render_aa myaa([4,2]); end






%% DEM w/ SCP, finding optimal K


% colormap info
cmin = 1;   cmax = 1.2;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');



% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

Wcot = makeOneRingWeights(mesh, 'dcp');

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
Wsave_dem_geoball = W;
W = Wsave_dem_geoball;  % use to restore if we computed earlier


W = Wsave_dem_geoball;  % use to restore if we computed earlier

% find optimal W scale value k  using SCP
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, W, x ), init_val/2, 2*init_val );
mesh.u = embedSCP(mesh, 'generalized', k*W );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, W, k, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
figure; plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('SCP w/ DEM Weights');


W = Wsave_dem_geoball;  % use to restore if we computed earlier
% scale weights by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    W(k,idx) = valence * W(k,idx);
end


% find optimal W scale value k  using SCP
init_val = 2;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, W, x ), init_val/2, 2*init_val );
mesh.u = embedSCP(mesh, 'generalized', k*W );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, W, k, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
figure; plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('SCP w/ valence-scaled DEM Weights');


W = Wsave_dem_geoball;  % use to restore if we computed earlier
%scale interior by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    if ( ~ isboundaryv(k) )
        W(k,idx) = valence * W(k,idx);
    end
end

% find optimal W scale value k  using SCP
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, W, x ), init_val/2, 2*init_val );
mesh.u = embedSCP(mesh, 'generalized', k*W );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, W, k, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
figure; plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('SCP w/ valence-scaled DEM Weights (interior)');





%% DEM w/ Conformal LLE, finding optimal K


% colormap info
cmin = 1;   cmax = 1.2;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');



% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);    F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

Wcot = makeOneRingWeights(mesh, 'dcp');

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
Wsave_dem_geoball = W;

W = Wsave_dem_geoball;  % use to restore if we computed earlier

% find optimal W scale value k using CLLE
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'lle', mesh, W, x, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh ), init_val/2, 2*init_val );
mesh.u = embedLLE(mesh.v, k*W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, W, k, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);


W = Wsave_dem_geoball;  % use to restore if we computed earlier
% scale weights by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    W(k,idx) = valence * W(k,idx);
end


% find optimal W scale value k  using SCP
init_val = 2;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'lle', mesh, W, x, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh ), init_val/2, 2*init_val );
mesh.u = embedLLE(mesh.v, k*W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, W, k, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ valence-scaled DEM Weights (QC:', num2str(sum(faceQC)-F),')']);


W = Wsave_dem_geoball;  % use to restore if we computed earlier
%scale interior by point valence
for k = 1:N
    idx = find(W(k,:));
    valence = numel(idx);
    if ( ~ isboundaryv(k) )
        W(k,idx) = valence * W(k,idx);
    end
end

% find optimal W scale value k  using SCP
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'lle', mesh, W, x, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh ), init_val/2, 2*init_val );
mesh.u = embedLLE(mesh.v, k*W, mesh.vidx(mesh.isboundaryv~=0), mesh.loops, mesh );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, W, k, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ valence-scaled DEM Weights (interior) (QC:', num2str(sum(faceQC)-F),')']);




%% failures on complex dog mesh


% colormap info
cmin = 1;   cmax = 1.5;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');
render_aa = 0;     % set to 1 to render hires AA images

mesh = readMesh('doghead_basehole.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);  F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh, options);
Wsave_dem_geoball = W;
W = Wsave_dem_geoball;  % use to restore if we computed earlier

% boundary constraints from [Mullen02]  (use real boundary)
mesh.u = embedLEM(mesh.v, W, mesh.vidx(mesh.isboundaryv~=0) );
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('LEM+RealBoundary - 2D');

mesh.u = embedSCP(mesh, 'generalized');
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('SCP');


k=1;
mesh.u = embedSCP(mesh, 'generalized', k*W);
figure;
plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);


faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ valence-scaled DEM Weights (interior) (QC:', num2str(sum(faceQC)-F),')']);


% find optimal W scale value k  using SCP
init_val = 10;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, x, W ), init_val/4, 4*init_val );
mesh.u = embedSCP(mesh, 'generalized', k*W );
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)),')']);


faceQC = faceDistortion(mesh, 'qc');
figure; plotMesh(mesh, 'UefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ valence-scaled DEM Weights (interior) (QC:', num2str(sum(faceQC)-F),')']);



%% energy landscape

mesh = readMesh('saddle.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');


% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
[W,err2,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
Wsave_dem_geoball = W;
W = Wsave_dem_geoball;  % use to restore if we computed earlier



X = [];
Cvals = [];    Dvals = [];   Avals = [];
kCvals = [];   kDvals = [];  kAvals = [];
Evals = [];
turnVals = [];
%for k = 1:0.025:40   % saddle, larger pass
for k = 1:0.002:20      % dogface
    mesh.u = embedSCP(mesh, 'generalized', k*W);
    turnVals = [turnVals, turningNumber(mesh.u(mesh.loops,:))];
    
    [c,d,a,t] = eConformal(mesh, k*W);
    E = FnOptLaplacianScale( 'scp', mesh, k, W );

    X = [X,k];
    kCvals = [kCvals,c];   kDvals = [kDvals,d];  kAvals = [kAvals,a];
    Evals = [Evals,E];
    
    plotMesh(mesh, 'Uefb');    
    drawnow; 
end

%turnVals2=floor(Evals/10);   % why doesn't turnVals == this ??
%plotvals = Evals-10*turnVals2;
plotvals = kCvals.^2;
figure;
hold all;
Cscale = mean(plotvals)/2;
plot(X,Cscale*turnVals);
plot(X,plotvals);
ylim([0,max(plotvals)*1.25])
hold off;
drawnow;


% sequence of plots
k=20;
k=14;
k=8.79;
k=6.92;
k=4.76;
k=3.63;
k=2.974;
mesh.u = embedSCP(mesh, 'generalized', k*W);
plotMesh(mesh, 'UefbO');
%figure; plotMesh(mesh, 'UefbO');
myaa;










