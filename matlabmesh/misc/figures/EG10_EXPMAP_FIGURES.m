%% fixed boundary param

% [RMS] this mesh shows LLE failure
mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;    N = size(P,1);
isboundaryv = mesh.isboundaryv;
%figure;   plotMesh(mesh,'efb');

boundaryUV = embedBoundary( mesh, 'circle' );




%% weight plots

interiorv=mesh.vidx(mesh.isboundaryv==0);
boundaryv=mesh.vidx(mesh.isboundaryv==1);
filterv = interiorv;


% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-4;
[W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
Wsave_dem_geoball = W;
W = Wsave_dem_geoball;  % use to restore if we computed earlier
optimization_init_val=10;

scrsz = get(0,'ScreenSize');
figure('Position',[200 scrsz(4)/2-200 scrsz(3)-200 scrsz(4)/3] );

subplot(1,3,1);
hist(nonzeros(W(:)),256);
title('weight histogram');

subplot(1,3,2);
hold on;
geodist = sqrt( Vu.^2 + Vv.^2 );
w_int = W(interiorv,:);
g_int = geodist(interiorv,:);
idx = find(w_int);
Wint = full(w_int(idx));
Gint = full(g_int(idx));
scatter(Gint,Wint);
w_bdr = W(boundaryv,:);
g_bdr = geodist(boundaryv,:);
idx = find(w_bdr);
Wbdr = full(w_bdr(idx));
Gbdr = full(g_bdr(idx));
scatter(Gbdr,Wbdr);
line([min(geodist(:)),max(geodist(:))],[0,0]);
title('weights vs geodesic distance');
hold off;

subplot(1,3,3);

allangles = [];
allweightsA = [];
for i = 1:N
    nbrs = find(W(i,:) > 0);
    K = numel(nbrs);
    u = full(Vu(i,nbrs));
    v = full(Vv(i,nbrs));
    angles = atan2(v,u);
    nbrangles = zeros(K,1);
    [angles,sidx] = sort(angles);
    for k = 1:K-1
        t = angles(k+1) - angles(k);
        nbrangles(sidx(k)) = nbrangles(sidx(k)) + t/2;        
        nbrangles(sidx(k+1)) = nbrangles(sidx(k+1)) + t/2;
    end
    t = angles(1)+2*pi - angles(K);
    nbrangles(sidx(K)) = nbrangles(sidx(K)) + t/2;        
    nbrangles(sidx(1)) = nbrangles(sidx(1)) + t/2;
    
    weights = full(W(i,nbrs));   
    allangles = [allangles; nbrangles];
    allweightsA = [allweightsA; weights'];
end

hold all;
line([min(allangles),max(allangles)],[0,0]);
scatter(allangles,allweightsA);
title('weights vs angular area');
hold off;

drawnow;






figure;
hold on;
hist(nonzeros(W(interiorv,:)),64);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','r')
hist(nonzeros(W(boundaryv,:)),64);
title('weight histogram');
hold off;




W = makeOneRingWeights(mesh,'dcp');


scrsz = get(0,'ScreenSize');
figure('Position',[200 scrsz(4)/2-200 scrsz(3)-200 scrsz(4)/3] );

subplot(1,3,1);
hist(nonzeros(W(:)),256);
title('weight histogram');

subplot(1,3,2);
hold on;
geodist = sqrt( Vu.^2 + Vv.^2 );
w_int = W(interiorv,:);
g_int = geodist(interiorv,:);
idx = find(w_int);
Wint = full(w_int(idx));
Gint = full(g_int(idx));
scatter(Gint,Wint);
w_bdr = W(boundaryv,:);
g_bdr = geodist(boundaryv,:);
idx = find(w_bdr);
Wbdr = full(w_bdr(idx));
Gbdr = full(g_bdr(idx));
scatter(Gbdr,Wbdr);
line([min(geodist(:)),max(geodist(:))],[0,0]);
title('weights vs geodesic distance');
hold off;

subplot(1,3,3);

allangles = [];
allweightsA = [];
for i = 1:N
    nbrs = find(W(i,:) > 0);
    K = numel(nbrs);
    u = full(Vu(i,nbrs));
    v = full(Vv(i,nbrs));
    angles = atan2(v,u);
    nbrangles = zeros(K,1);
    [angles,sidx] = sort(angles);
    for k = 1:K-1
        t = angles(k+1) - angles(k);
        nbrangles(sidx(k)) = nbrangles(sidx(k)) + t/2;        
        nbrangles(sidx(k+1)) = nbrangles(sidx(k+1)) + t/2;
    end
    t = angles(1)+2*pi - angles(K);
    nbrangles(sidx(K)) = nbrangles(sidx(K)) + t/2;        
    nbrangles(sidx(1)) = nbrangles(sidx(1)) + t/2;
    
    weights = full(W(i,nbrs));   
    allangles = [allangles; nbrangles];
    allweightsA = [allweightsA; weights'];
end

hold all;
line([min(allangles),max(allangles)],[0,0]);
scatter(allangles,allweightsA);
title('weights vs angular area');
hold off;

drawnow;




%% weight histograms

mesh = readMesh('bunny.obj');


interiorv=mesh.vidx(mesh.isboundaryv==0);
boundaryv=mesh.vidx(mesh.isboundaryv==1);
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
%options.nbrtype = 'geoball';
%options.nbrsize = 2.1*edgelen_avg;
options.nbrtype = 'k';
options.nbrsize = 20;

xlo = -0.5;     xhi = 0.8;
ylo = 0;        yhi = 800;
show_boundary = 0;
show_interior = 1;
BINS = xlo:0.025:xhi;

W = makeOneRingWeights(mesh,'dcp',1);

scrsz = get(0,'ScreenSize');
figure('Position',[900 scrsz(4)/2-200 700 700] );
set(gcf,'DefaultAxesFontSize',20);
hold on;
if show_interior
    hist(nonzeros(W(interiorv,:)),BINS);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','r')
end
if show_boundary
    hist(nonzeros(W(boundaryv,:)),BINS);
end
ylim([ylo,yhi]); xlim([xlo,xhi]);
title('Wcotan');
hold off;

%showvals = [10e-1,10e-3,10e-5];
showvals = [10e-5];
for k = 1:numel(showvals)
    options.regtol = showvals(k);
    [W,err,Vu,Vv] = makeExpMapWeights(mesh, options);
    scrsz = get(0,'ScreenSize');
    figure('Position',[900 scrsz(4)/2-200 700 700] );
    set(gcf,'DefaultAxesFontSize',20);    hold on;
    if show_interior
        hist(nonzeros(W(interiorv,:)),BINS);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','EdgeColor','r')
    end
    if show_boundary
        hist(nonzeros(W(boundaryv,:)),BINS);
    end
    ylim([ylo,yhi]); xlim([xlo,xhi]);
    title(['Wtangent, eps:',num2str(options.regtol),' recons err:',num2str(sum(err))]);
    hold off;
end




%% DEM speed

mesh = readMesh('dogface.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;

mesh = readMesh('doghead_basehole.obj');
P = mesh.v;  Pn = mesh.n;   N = size(P,1);
isboundaryv = mesh.isboundaryv;

% |expmapnbrs| ~ 15
for k = 1:5
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = 0.97*edgelen_avg;
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
end


% |expmapnbrs| ~ 30
for k = 1:5
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = 1.85*edgelen_avg;
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
end


% |expmapnbrs| ~ 60
for k = 1:5
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = 3*edgelen_avg;
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
end



% |Ni| = 15
for k = 1:5
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'k';
    options.nbrsize = 15;
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
end


% |Ni| = 30
for k = 1:5
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'k';
    options.nbrsize = 30;
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
end


% |Ni| = 60
for k = 1:5
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'k';
    options.nbrsize = 60;
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);
end


useK = [15,30,60]
for i=1:numel(useK)
    % DEM weights, G=2.1*edgelenavg
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    options.weightmode = 'optimal2';
    options.nbrtype = 'k';
    options.nbrsize = useK(i);
    [W,err,Vu,Vv] = makeExpMapWeights(mesh.v, mesh.n, options);

    % find optimal W scale value k  using SCP
    [k,fval,exitflag] = fminbnd( @(x) FnOptLaplacianScale( 'scp', mesh, x, W ), 5, 50 );
    mesh.u = embedSCP(mesh, 'generalized', k*W );
    fprintf('Optimal k: %f \n', k);
    figure; plotMesh(mesh, 'UefbO');

    % average time to solve (over 5 solves)
    tic;
    for t = 1:5
        mesh.u = embedSCP(mesh, 'generalized', k*W );
    end
    solve5 = toc;
    fprintf('solve time (%d nbrs): %f\n', useK(i), solve5/5);
end



%% euclidean vs tangent-space weights

mesh = readMesh('4bump.obj');
mesh = readMesh('bump_ref.obj');
mesh = readMesh('bump_mc.obj');

% DEM weights, G=2.1*edgelenavg
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.1*edgelen_avg;
options.regtol = 10e-3;
Wtan = makeExpMapWeights(mesh, options);

% DEM weights, G=2.1*edgelenavg
options.weightmode = 'optimal3';
W3d = makeExpMapWeights(mesh, options);


% LLE
mesh.u = embedLLE(mesh.v, Wtan );
figure; plotMesh(mesh, 'Uefb'); title('LLE - Wtan');
plotMeshGL(mesh);
pause(10);   % give plotMeshGL time to read file... (fix this??)

mesh.u = embedLLE(mesh.v, W3d );
figure; plotMesh(mesh, 'Uefb'); title('LLE - W3d');
plotMeshGL(mesh);



% LEM
mesh.u = embedLEM(mesh.v, Wtan, mesh.loops{1} );
figure; plotMesh(mesh, 'Uefb'); title('LLE - Wtan');
plotMeshGL(mesh);
pause(10);   % give plotMeshGL time to read file... (fix this??)

mesh.u = embedLEM(mesh.v, W3d, mesh.loops{1} );
figure; plotMesh(mesh, 'Uefb'); title('LLE - W3d');
plotMeshGL(mesh);


%% reconstruction error plots

mesh = readMesh('face_10k_holes.pobj');
plotMesh(mesh,'efb');
mesh.u = embedSCP(mesh);
plotMesh(mesh,'uefb');

N = numel(mesh.vidx);
mesh.v = [mesh.u,zeros(N,1)];

Wcot = makeOneRingWeights(mesh,'dcp');
Wcot2 = Wcot;
Wcot2(1:N+1:N*N) = -sum(Wcot2,2);
Rerr = (Wcot2*mesh.v(:,1)).^2 + (Wcot2*mesh.v(:,2)).^2;

tmp = mesh;
tmp.v = [mesh.u,50*Rerr];
tmp.n = estimateNormal(tmp);
plotMesh(tmp,'efbl');
view(45,45)


pointset = makePointSet(mesh.v,mesh.n,mesh.u);
edgelen_avg = edgeStats(pointset.v, pointset.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
options.regtol = 10e-3;
Wtan = makeExpMapWeights(pointset, options);


Wtan2 = Wtan;
Wtan2(1:N+1:N*N) = -sum(Wtan2,2);
Rerr = (Wtan2*mesh.v(:,1)).^2 + (Wtan2*mesh.v(:,2)).^2;



tmp = mesh;
tmp.v = [mesh.u,50*Rerr];
tmp.n = estimateNormal(tmp);
plotMesh(tmp,'efbl');
view(45,45)




