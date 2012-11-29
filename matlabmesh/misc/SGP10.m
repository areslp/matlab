

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis of eigenspectrum of laplace-beltrami operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mesh = readMesh('doghead.obj');
mesh = readMesh('sphere.obj');
mesh = readMesh('sphere_semireg_1000.obj');
mesh = readMesh('sphere_semireg_2000.obj');
mesh = readMesh('sphere_semireg_4000.obj');
mesh = readMesh('sphere_graphite_300.obj','n');
mesh = readMesh('sphere_graphite_2000.obj','n');

N = size(mesh.v,1);
maxE = 100;
figure; plotMesh(mesh,'efb');
eigsopts = struct('disp',0);

% COTAN

Wcot = makeOneRingWeights(mesh, 'dcp');
Wcot = -Wcot;
Wcot(1:N+1:N*N) = -sum(Wcot,2);  % set W(i,i) = -sum_j W(i,j)
Amix = sparse(N,N);
Amix(1:N+1:N*N) = vertexArea(mesh, [], 'mixed');
dNoArea = real(eigs(Wcot, maxE, -1e-5,eigsopts));        if dNoArea(1)>dNoArea(end) dNoArea = dNoArea(maxE:-1:1); end;
figure('Name','Cotan'); plot(dNoArea(1:maxE));
d = real(eigs(Wcot, Amix, maxE, -1e-5,eigsopts));        if d(1)>d(end) d = d(maxE:-1:1); end;
figure('Name','Area-Weighted Cotan'); plot(d(1:maxE));

% NORMALIZED COTAN

WcotN = makeOneRingWeights(mesh, 'dcp', 1);
WcotN = -WcotN;
WcotN(1:N+1:N*N) = -sum(WcotN,2);  % set W(i,i) = -sum_j W(i,j)
dN = real(eigs(WcotN, maxE, -1e-5,eigsopts));
figure('Name','Normalized Cotan'); plot(dN(1:maxE));
dN = real(eigs(WcotN, Amix, maxE, -1e-5,eigsopts));
figure('Name','Area-Weighted Normalized Cotan'); plot(dN);

% UNIFORM
Wuni = makeOneRingWeights(mesh, 'uniform');
Wuni = -Wuni;
Wuni(1:N+1:N*N) = -sum(Wuni,2);  % set W(i,i) = -sum_j W(i,j)
dUni = real(eigs(Wuni, maxE, -1e-5,eigsopts));
if dUni(1)>dUni(end) dUni = dUni(maxE:-1:1); end;
figure('Name','Uniform'); plot(dUni(1:maxE));
dUni = real(eigs(Wuni, Amix, maxE, -1e-5,eigsopts));
if dUni(1)>dUni(end) dUni = dUni(maxE:-1:1); end;
figure('Name','Area-Weighted Uniform'); plot(dUni(1:maxE));


% one-ring optimal2

options = struct('regtol',10e-1);
Wopt2 = makeOneRingWeights(mesh, 'optimal2', 0, 'uniform',options);
Wopt2 = -Wopt2;
Wopt2(1:N+1:N*N) = -sum(Wopt2,2);  % set W(i,i) = -sum_j W(i,j)
dOpt2 = real(eigs(Wopt2, maxE, -1e-5,eigsopts)); 
figure('Name','One-Ring Optimal2'); plot(dOpt2);
dOpt2 = real(eigs(Wopt2, Amix, maxE, -1e-5,eigsopts)); 
figure('Name','Area-Weighted One-Ring Optimal2'); plot(dOpt2);
  

% GEOBALL HEATKERNEL

edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'gaussian3';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
%options.heath = 0.0025;
options.heath = 0.0141;
[Wheat3,err,vU,vV] = makeExpMapWeights(mesh, options);
Wheat3 = -Wheat3;
Wheat3(1:N+1:N*N) = -sum(Wheat3,2);  % set W(i,i) = -sum_j W(i,j)
dH3 = real(eigs(Wheat3, maxE, -1e-5,eigsopts)); 
figure('Name',['Geoball HeatK, heath =',num2str(options.heath)]); plot(dH3);

ptareaopt = struct('u',vU,'v',vV,'radius', sqrt(options.heath) );
[tanA,isb] = pointArea(mesh.v, [], 'uvDelArea', ptareaopt);
invD = sparse(1:N,1:N, tanA);
dH3A = real(eigs(Wheat3, invD, maxE, -1e-5,eigsopts));
figure('Name','delaunay area-weighted HeatK'); plot(dH3A);


% GEOBALL EXPMAP

edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
%options.regtol=10e-5;
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.0*edgelen_avg;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
fprintf('recons err is %f\n',sum(err));
Wtan = -Wtan;
Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)
d2 = real(eigs(Wtan, maxE, -1e-5,eigsopts));
figure('Name','Geoball WTan'); plot(d2);

ptareaopt = struct('u',vU,'v',vV);
[tanA,isb] = pointArea(mesh.v, [], 'uvVoronoi', ptareaopt);
invD = sparse(1:N,1:N, tanA);
d2 = real(eigs(Wtan, invD, maxE, -1e-5,eigsopts));
figure('Name','delaunay area-weighted geoball WTan'); plot(d2);


% K-NBR EXPMAP

edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 20;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
Wtan = -Wtan;
Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)
d2 = real(eigs(Wtan, maxE, -1e-5,eigsopts));
figure('Name','K=20 WTan'); plot(d2);

ptareaopt = struct('u',vU,'v',vV);
[tanA,isb] = pointArea(mesh.v, [], 'uvVoronoi', ptareaopt);
invD = sparse(1:N,1:N, tanA);
d2 = real(eigs(Wtan, invD, maxE, -1e-5,eigsopts));
figure('Name','delaunay area-weighted K=20 WTan'); plot(d2);



d2 = eig(full(Wtan'*Wtan));
figure('Name','Wtan*Wtan'); plot(real(d2(1:maxE)));
dN = eig(full(WcotN'*WcotN));
figure('Name','Normalized Wcot*Wcot'); plot(dN(1:maxE));
d = eig(full(Wcot'*Wcot),full(Amix));
figure('Name','Area-Weighted Wcot*Wcot'); plot(d(1:maxE));
d = eig(full(Wcot'*Wcot);
figure('Name','Wcot*Wcot'); plot(d(1:maxE));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% overlay DEM param with DCP param (circle boundary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tmp = readMesh('patch.obj');
tmp = clipEars(tmp);
boundaryUV = embedBoundary(tmp, 'circle');

options = struct('regtol',10^-3);
useW = makeOneRingWeights(tmp, 'optimal2', 0, 'uniform', options);
tmp.u = embedInterior(tmp,boundaryUV,useW);
figure;plotMesh(tmp,'uefb');


options = struct('regtol',10^-2.2);
useW = makeOneRingWeights(tmp, 'optimal3', 0, 'uniform', options);
tmp.u = embedInterior(tmp,boundaryUV,useW);
figure;plotMesh(tmp,'uefb');


figure; hold all;
tmp.u = embedInterior(tmp,boundaryUV,useW);
plotMesh(tmp,'ue');
tmp.u = embedInterior(tmp,boundaryUV,makeOneRingWeights(tmp,'dcp'));
plotMesh(tmp,'ue');
hold off;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (more) Analysis of eigenspectrum of laplace-beltrami operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
fprintf('recons err is %f\n',sum(err));
ptareaopt = struct('u',vU,'v',vV);
[tanA,isb] = pointArea(mesh.v, [], 'uvVoronoi', ptareaopt);

Wout = sparse(N,N);
for k = 1:N
   G = -Wtan(k,:); 
   G(k) = -sum(G .* tanA') / tanA(k) + G(k);
   Wout(k,:) = G;
end
invD = sparse(1:N,1:N, 1./tanA);
d2 = real(eigs(Wout, invD, maxE, -1e-5,eigsopts));
figure('Name','gnl area-weighted WTan'); plot(d2);

WtanL = -Wtan;
WtanL(1:N+1:N*N) = -sum(WtanL,2);  % set W(i,i) = -sum_j W(i,j)
d2 = real(eigs(WtanL, maxE, -1e-5,eigsopts));
figure('Name','WTan'); plot(d2);




edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
options.weightmode = 'optimal2b';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
fprintf('recons err is %f\n',sum(err));
WtanL = -Wtan;
WtanL(1:N+1:N*N) = -sum(WtanL,2);  % set W(i,i) = -sum_j W(i,j)
d2 = real(eigs(WtanL, maxE, -1e-5,eigsopts));
figure('Name','WTan'); plot(d2);



edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
ptareaopt = struct('u',vU,'v',vV);
[tanA,isb] = pointArea(mesh.v, [], 'uvVoronoi', ptareaopt);
options.areas = tanA;
options.weightmode = 'optimal2b';
[Wtan2b,err,vU,vV] = makeExpMapWeights(mesh, options);

WtanL = -Wtan2b;
WtanL(1:N+1:N*N) = -sum(WtanL,2);  % set W(i,i) = -sum_j W(i,j)
invD = sparse(1:N,1:N, tanA);
d2 = real(eigs(WtanL, invD, maxE, -1e-5,eigsopts));
figure('Name','optimal2b WTan'); plot(d2);









check_vals = [2.0,3.0,4.0];
evals = [];
avals = [];
for k = 1:numel(check_vals) 
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = check_vals(k)*edgelen_avg;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    ptareaopt = struct('u',vU,'v',vV);
    [tanA,isb] = pointArea(mesh.v, [], 'uvDelArea', ptareaopt);
    fprintf('recons err for nbrsize %f is %f\n',options.nbrsize,sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)
    d2 = real(eigs(Wtan, maxE, -1e-5,eigsopts));
    evals(k,:) = d2;
    avals(k,:) = tanA;
end

figure; hold all;
for k = 1:numel(check_vals) 
    r = check_vals(k)*edgelen_avg;
    E = evals(k,:);
    A = avals(k,:);
%    plot(E / (pi*r^2));
    plot(E / mean(A));
end
hold off; drawnow;


meshes = {'sphere_semireg_1000.obj','sphere_semireg_2000.obj','sphere_semireg_4000.obj'};
evals = [];
rvals = [];
for k = 1:numel(meshes)
    mesh = readMesh(meshes{k});
    N = size(mesh.v,1);
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    rvals(k) = edgelen_avg;
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = 3.0*edgelen_avg;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    fprintf('recons err for nbrsize %f is %f\n',options.nbrsize,sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)
    d2 = real(eigs(Wtan, maxE, -1e-5,eigsopts));
    evals(k,:) = d2;
end


figure; hold all;
for k = 1:numel(check_vals) 
    r = 3.0*rvals(k);
    plot(evals(k,:) / (pi*r^2));
end
hold off; drawnow;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare error in DEM vs analytic on hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE you need to mex util\dijkstra*.cpp for this to work

writedata = 1;
writedata = 0;

% construct hemisphere from sphere submesh
mesh = readMesh('sphere_semireg_2000.obj');
mesh.v = normalize(mesh.v);
radius = 1;

tzs = [ mesh.v(mesh.f(:,1),3), mesh.v(mesh.f(:,2),3), mesh.v(mesh.f(:,3),3) ];
hemifaces = find(tzs(:,1)>0 & tzs(:,2)>0 & tzs(:,3)>0);
hemimesh = subMesh(mesh, hemifaces);
figure;plotMesh(hemimesh);
northpole = [0,0,1];

% read hemisphere
hemimesh = readMesh('hemisphere_quadreg.obj');
hemimesh = readMesh('hemisphere_graphite_1000.obj', 'n');
hemimesh = readMesh('hemisphere_graphite_2000.obj', 'n');
hemimesh = readMesh('hemisphere_graphite_4000.obj', 'n');

hemimesh.v = normalize(hemimesh.v);
radius = 1;
figure;plotMesh(hemimesh);
northpole = [0,0,1];

[tan1,tan2] = tangentFrame(normalize(northpole));

meshuv = zeros(numel(hemimesh.vidx),2);
for vi = 1:numel(hemimesh.vidx)
    meshuv(vi,:) = sphereNormCoords(radius, northpole, hemimesh.v(vi,:), tan1,tan2);
end
hemimesh.u = meshuv;
meshG = vmag(meshuv);   meshT = atan2(meshuv(:,2),meshuv(:,1));
if (writedata) csvwrite('meshuv.csv', [meshG,meshT,meshuv]); end
[tmp,idx] = sort(meshG);



hemips = [];
hemips.v = [ northpole; hemimesh.v ];
hemips.n = normalize(hemips.v);

options = [];
options.graphnbrs = 8;
options.stopcrit = 'none';
options.stopparam = inf;
options.upwindavg = 0;
expmapuv = expmap(hemips.v, hemips.n, 1, options);
expmapuv = expmapuv(2:end,:);
rotangle = meshT(1)-atan2(expmapuv(1,2),expmapuv(1,1));
rotmat = [cos(rotangle),sin(rotangle); -sin(rotangle),cos(rotangle)];
expmapuv = mvmul(expmapuv,rotmat)';
expmapG = vmag(expmapuv);  expmapT = atan2(expmapuv(:,2),expmapuv(:,1));
if (writedata) csvwrite('expmapuv.csv', [expmapG,expmapT,expmapuv]); end

options.upwindavg = 1;
upwinduv = expmap(hemips.v, hemips.n, 1, options);
upwinduv = upwinduv(2:end,:);
rotangle = meshT(1)-atan2(upwinduv(1,2),upwinduv(1,1));
rotmat = [cos(rotangle),sin(rotangle); -sin(rotangle),cos(rotangle)];
upwinduv = mvmul(upwinduv,rotmat)';
upwindG = vmag(upwinduv);  upwindT = atan2(upwinduv(:,2),upwinduv(:,1));
if (writedata) csvwrite('upwinduv.csv', [upwindG,upwindT,upwinduv]); end


figure('Name','Dist Error');hold all;
expmapGerr = abs(meshG-expmapG);
plot(expmapGerr(idx));      
upwindGerr = abs(meshG-upwindG);
plot(upwindGerr(idx));      
hold off; drawnow;

figure('Name','Angle Error');hold all;
expmapTerr=angledist(meshT,expmapT);
plot(abs(expmapTerr(idx)));
upwindTerr=angledist(meshT,upwindT);
plot(abs(upwindTerr(idx)));
hold off; drawnow;


figure('Name','L2 Error');hold all;
expmapUVerr=vmag2(meshuv-expmapuv);
plot(abs(expmapUVerr(idx)));
upwindUVerr=vmag2(meshuv-upwinduv);
plot(abs(upwindUVerr(idx)));
hold off; drawnow;





figure('Name','Dist Error');hold all;
expmapGerr = abs(meshG-expmapG);
plot(sort(expmapGerr));      
upwindGerr = abs(meshG-upwindG);
plot(sort(upwindGerr));      
hold off; drawnow;

figure('Name','Angle Error');hold all;
expmapTerr=angledist(meshT,expmapT);
plot(abs(sort(expmapTerr)));
upwindTerr=angledist(meshT,upwindT);
plot(abs(sort(upwindTerr)));
hold off; drawnow;

figure('Name','L2 Error');hold all;
expmapUVerr=vmag2(meshuv-expmapuv);
plot(abs(sort(expmapUVerr)));
upwindUVerr=vmag2(meshuv-upwinduv);
plot(abs(sort(upwindUVerr)));
hold off; drawnow;



figure; newplot; hold all; 
hemimesh.u = meshuv;
plotMesh(hemimesh,'ue');
hemimesh.u = expmapuv;
plotMesh(hemimesh,'ue');
hold off; drawnow;


figure; newplot; hold all; 
hemimesh.u = meshuv;
plotMesh(hemimesh,'ue');
hemimesh.u = upwinduv;
plotMesh(hemimesh,'ue');
hold off; drawnow;


options = [];
options.weightmode = 'uniform';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
%options.heath = 0.0025;
options.heath = 0.0141;
[Wheat3,err,vU,vV] = makeExpMapWeights(mesh, options);
Wheat3 = -Wheat3;
Wheat3(1:N+1:N*N) = -sum(Wheat3,2);  % set W(i,i) = -sum_j W(i,j)
dH3 = real(eigs(Wheat3, maxE, -1e-5,eigsopts)); 
figure('Name',['Geoball HeatK, heath =',num2str(options.heath)]); plot(dH3);





