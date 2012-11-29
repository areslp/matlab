


%% optDNCP - saddle

% colormap info
cmin = 1;   cmax = 1.1;     clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');
render_aa = 0;     % set to 1 to render hires AA images

mesh = readMesh('saddle.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

W = makeOneRingWeights(mesh, 'dcp');
optimization_init_val = 1;



mesh.u = embedSCP(mesh, 'generalized');
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;


constrain_vtx1 = 1;
constrain_vtx2 = 102;

mesh.u = embedDNCP(mesh, W, constrain_vtx1, constrain_vtx2);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;


% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;




constrain_vtx1 = 1;
constrain_vtx2 = 4;

mesh.u = embedDNCP(mesh, W, constrain_vtx1, constrain_vtx2);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;


% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;



%% optDNCP - other examples


cmin = 1;   cmax = 1.3;

mesh = readMesh('dogface.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;


W = makeOneRingWeights(mesh, 'dcp');
optimization_init_val = 1;


constrain_vtx1 = 1;
constrain_vtx2 = 320;

mesh.u = embedDNCP(mesh, W, constrain_vtx1, constrain_vtx2);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;


% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;






cmin = 1;   cmax = 1.7;

mesh = readMesh('doghead_basehole.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;


W = makeOneRingWeights(mesh, 'dcp');
optimization_init_val = 1;

constrain_vtx1 = 1118;
constrain_vtx2 = 1545;

mesh.u = embedDNCP(mesh, W, constrain_vtx1, constrain_vtx2);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;


% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
consv = mesh.u([constrain_vtx1,constrain_vtx2],:);
scatter(consv(:,1), consv(:,2), 100, 'ms', 'filled');
hold off; drawnow;




%% optDNCP - different weights


W = makeOneRingWeights(mesh, 'uniform');
optimization_init_val = 1;
constrain_vtx1=-1;
constrain_vtx2=-1;
cmin = 1;   cmax = 1.2; 

% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO');
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;


W = makeOneRingWeights(mesh, 'uniform', 1);
optimization_init_val = 4;

% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;



W = makeOneRingWeights(mesh, 'dcp', 1);
optimization_init_val = 4;

% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'UefbO');
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;



% meanvalue
W = makeOneRingWeights(mesh, 'meanvalue');
cmin = 1;   cmax = 1.1; 


% find optimal W scale value k  using DNCP
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W), 0.001, 10 );
mesh.u = embedDNCP(mesh, k*W );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
%plotMesh(mesh, 'uefbO', faceQC);
plotMesh(mesh, 'uefbO');
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;




% dap
W = makeOneRingWeights(mesh, 'dap');
cmin = 1;   cmax = 1.1; 


% find optimal W scale value k  using DNCP
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W), 0.001, 10 );
mesh.u = embedDNCP(mesh, k*W );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
%plotMesh(mesh, 'uefbO', faceQC);
plotMesh(mesh, 'uefbO');
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;




mesh = readMesh('dogface.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');



W = makeOneRingWeights(mesh, 'uniform');
optimization_init_val = 1;
constrain_vtx1=-1;
constrain_vtx2=-1;
cmin = 1;   cmax = 1.2; 

% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO');
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;









mesh = readMesh('bump_ref.obj');
mesh.n = estimateNormal(mesh);
P = mesh.v;  Pn = mesh.n;   N = size(P,1);   F = size(mesh.f,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');



W = makeOneRingWeights(mesh, 'uniform');
optimization_init_val = 1;
constrain_vtx1=-1;
constrain_vtx2=-1;
cmin = 1;   cmax = 1.2; 

% find optimal W scale value k  using DNCP
init_val = optimization_init_val;
options = [];
options.TolX = 0.0000001;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'dncp', mesh, x, W, constrain_vtx1, constrain_vtx2 ), init_val/2, 2*init_val );
mesh.u = embedDNCP(mesh, k*W, constrain_vtx1, constrain_vtx2 );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*W );
E = FnOptLaplacianScale( 'dncp', mesh, k, W, constrain_vtx1, constrain_vtx2 );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
fprintf('Optimal k: %f \n', k);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO');
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;






%% robustness

clo = [.3,.3,1];    chi = [0,1,0];
load('myjet.mat');


mesh = readMesh('bunnyhead1.obj');
mesh.v = [mesh.v(:,1), mesh.v(:,2), mesh.v(:,3)];
P = mesh.v;  Pn = mesh.n;    N = size(P,1);
isboundaryv = mesh.isboundaryv;
figure;   plotMesh(mesh,'efb');

Wcot = makeOneRingWeights(mesh, 'dcp');
Wuni = makeOneRingWeights(mesh, 'uniform');
optimization_init_val = 1;
constrain_vtx1 = -1;
constrain_vtx2 = -1;

cmin = 1;   cmax = 4;    
Wn = Wcot;  Wn(Wn(:)>0) = 0;
Wnsum = abs(sum(Wn,2));
plotMesh(mesh,'efb',Wnsum);
%colormap(myjet); caxis([cmin,cmax]);
colormap(summer);  caxis([cmin,cmax]);

cmin = 1;   cmax = 1.2;    
mesh.u = embedSCP(mesh, 'generalized', Wcot);
faceQC = faceDistortion(mesh, 'qc');
figure; hold all;
plotMesh(mesh, 'uefbO', faceQC);
colormap(myjet); caxis([cmin,cmax]);
title(['C-LLE w/ DEM Weights (QC:', num2str(sum(faceQC)-F),')']);
hold off; drawnow;


% construct mixed cotan/uniform weight matrix using
% some 3D balls around specific vertices
Wmix = Wcot;
Vchanged = zeros(N,1);
for k = 1:N
    dist1 = vmag2(mesh.v(k,:) - mesh.v(15,:));
    dist2 = vmag2(mesh.v(k,:) - mesh.v(69,:));
    if min(dist1,dist2) < 0.16
        Wmix(k,:) = Wuni(k,:);
        Vchanged(k) = 1;
    end
end
figure;
plotMesh(mesh,'efb',Vchanged);
colormap('summer');

% find optimal W scale value k  using SCP
init_val = 1;
[k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, x, Wmix ), init_val/2, 2*init_val );
mesh.u = embedSCP(mesh, 'generalized', k*Wmix );
[ EConformal, EDirichlet, Area, Turning ] = eConformal( mesh, k*Wmix );
E = FnOptLaplacianScale( 'scp', mesh, k, Wmix );
fprintf('k: %f   E: %f   C: %f   turning: %f\n', k, E, EConformal, Turning);
figure; plotMesh(mesh, 'UefbO', faceDistortion(mesh, 'qc'));
colormap(myjet); caxis([cmin,cmax]);
title('SCP w/ valence-scaled DEM Weights');