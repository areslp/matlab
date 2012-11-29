
%% L-B operator tests on plane

% fixed seed PRNG so we can repeat computations
rand('twister',1702925028);

% generate points on plane
N = 500;
%uv = genpoints(N,2,'grid_unit_square');
uv = genpoints(N,2,'stratified_unit_square');
%uv = genpoints(N,2,'uniform_unit_square');
uv = vadd(2*uv,[-1,-1]);
N = size(uv,1);
TRI = delaunay(uv(:,1), uv(:,2));
mesh = makeMesh( [uv,zeros(N,1)], TRI, repmat([0,0,1],N,1) );
figure;
plotMesh(mesh);
P = mesh.v;
Pn = mesh.n;
n = size(P,1);
%figure;
%scatter3(P(:,1),P(:,2),P(:,3));
% only want to use real interior verts, because estimation near boundary is
%  very bad in most cases...
dthresh=0.6;
interiorv = (mesh.v(:,1) > -dthresh) &  (mesh.v(:,1) < dthresh)  &   (mesh.v(:,2) > -dthresh)  &    (mesh.v(:,2) < dthresh);

% L-B functions on planes

f = 2*P(:,1);           %  f = 2x
LBf = zeros(size(f));   %  f_xx + f_yy = 0


f = P(:,1).^2;          %  f = x^2
LBf = 2*ones(size(f));  %  f_xx + f_yy = 2

f = exp(P(:,1)+P(:,2));  % f = e^(x+y)
LBf = 2*f;               % f_xx + f_yy = 2e^(x+y)





for nbrsize=6:4:30

% DEM weights w/o using mesh
emoptions = [];
emoptions.weightmode = 'optimal2';
%emoptions.nbrtype = 'geoball';
%emoptions.nbrsize = 0.2;
%emoptions.nbrsize = 0.025*nbrsize;
emoptions.nbrtype = 'k';
%emoptions.nbrsize = 10;
emoptions.nbrsize = nbrsize;
emoptions.silent = 1;
[W,err2,Vu,Vv] = makeExpMapWeights(P, Pn, emoptions);

% estimate radius of integration disc  (can directly use nbrsize for geoball...)
rEst = zeros(n,1);
for k = 1:n
    rEst(k) = max(vmag2([full(nonzeros(Vu(k,:))),full(nonzeros(Vv(k,:)))]));
end
intR = mean(rEst);
meanW = mean(nonzeros(W(:)));
    
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,isboundary] = pointArea(P, [], 'uvVoronoi', aoptions);
Adelsum = pointArea(P, [], 'uvDelArea', aoptions);

WA = W;     EM_A = A;
for k = 1:n
    WA(k,:) = W(k,:) .* A';
    EM_A(k) = (abs(W(k,:))>0)*A + A(k);
end
D = sparse(1:n,1:n,sum(W,2),n,n);
DA = sparse(1:n,1:n,sum(WA,2),n,n);

%LBl = -((DA-WA)*f)/2  ./ ( sum(WA,2) ) ./ EM_A / 0.1;
LBl = -((DA-WA)*f)/2   ./ ( sum(WA,2) ) ./ EM_A / meanW;
%LBl = -((DA-WA)*f)/2  ./ ( pi*realR^2*meanW ) ./ (pi*realR^2) / (meanW*2);
%LBl = -((DA-WA)*f)/2  ./ (EM_A * meanW) ./ EM_A / (meanW * 2);
%LBl = -((DA-WA)*f)/2   * 10000 * 2.4; % ./ (4*pi*0.2^2);
%LBl = -((DA-WA)*f)/2   / (mean(EM_A)*mean(A)); % ./ (4*pi*0.2^2);

%LBl = -(D-W)*f ./ (A) * 2;
%LBl = -(D-W)*f ./ (Adelsum.^2) * 2;
%LBl = -(D-W)*f ./ (EM_A.^2) * 2;

if norm(LBf) == 0
    L2 = norm(LBf-LBl);
    L2int = norm(LBf(interiorv)-LBl(interiorv));
    Lmax = max(abs(LBf - LBl));
    Lmaxint = max(abs(LBf(interiorv) - LBl(interiorv)));
    fprintf('#P=%d,sz=%2.2f    L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', N, emoptions.nbrsize, L2,L2int, Lmax,Lmaxint);
else 
    L2 = norm(LBf - LBl) / norm(LBf);
    L2int = norm(LBf(interiorv)-LBl(interiorv)) / norm(LBf(interiorv));
    Lmax = max(abs(LBf - LBl)) / max(abs(LBf));
    Lmaxint = max(abs(LBf(interiorv) - LBl(interiorv))) / max(abs(LBf(interiorv)));
    fprintf('#P=%d,sz=%2.2f    Normalized L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f  ratio: %f\n', N, emoptions.nbrsize, L2,L2int, Lmax,Lmaxint, mean(LBf(interiorv)./LBl(interiorv)));
%    fprintf('          meanW: %f   estR: %f  \n', meanW, intR );
end

end

subplot(1,2,1);
tmpv = [mesh.v(:,1),mesh.v(:,2),abs(LBl-LBf)];
tmpmesh = makeMesh(tmpv, mesh.f, mesh.n);
plotMesh(tmpmesh);

surf(reshape(abs(LBl-LBf), sqrt(N), sqrt(N)));



%
% cotan L-B
%
Wcot = makeOneRingWeights(mesh, 'dcp');
Dcot = sparse(1:n,1:n,sum(Wcot,2),n,n); 
Amix = vertexArea(mesh,[],'mixed');
LBlcot = -(Dcot-Wcot)*f ./ (2*Amix);

if norm(LBf) == 0
    L2 = norm(LBf-LBlcot);
    L2int = norm(LBf(interiorv)-LBlcot(interiorv));
    Lmax = max(abs(LBf - LBlcot));
    Lmaxint = max(abs(LBf(interiorv) - LBlcot(interiorv)));
    fprintf('L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', L2,L2int, Lmax,Lmaxint);
else 
    L2 = norm(LBf - LBlcot) / norm(LBf);
    L2int = norm(LBf(interiorv)-LBlcot(interiorv)) / norm(LBf(interiorv));
    Lmax = max(abs(LBf - LBlcot)) / max(abs(LBf));
    Lmaxint = max(abs(LBf(interiorv) - LBlcot(interiorv))) / max(abs(LBf(interiorv)));
    fprintf('#P=%d,COT Normalized L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', N, L2,L2int, Lmax,Lmaxint);
end




%
% expmap heat-kernel L-B
%
% DEM weights w/o using mesh
emoptions = [];
emoptions.weightmode = 'gaussian';
emoptions.nbrtype = 'geoball';
emoptions.nbrsize = 0.3;
emoptions.heath = 0.04/4;
[Weg,err2,Vu,Vv] = makeExpMapWeights(P, Pn, emoptions);
Aeg = vertexArea(mesh,[],'onering') / 3;
WAeg = Weg;
for k = 1:n
    WAeg(k,:) = Weg(k,:) .* Aeg';
end
Deg = sparse(1:n,1:n,sum(WAeg,2),n,n);
LBleg = -(Deg-WAeg)*f  / (4*pi*emoptions.heath^2);

if norm(LBf) == 0
    L2 = norm(LBf-LBleg);
    L2int = norm(LBf(interiorv)-LBleg(interiorv));
    Lmax = max(abs(LBf - LBleg));
    Lmaxint = max(abs(LBf(interiorv) - LBleg(interiorv)));
    fprintf('L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', L2,L2int, Lmax,Lmaxint);
else 
    L2 = norm(LBf - LBleg) / norm(LBf);
    L2int = norm(LBf(interiorv)-LBleg(interiorv)) / norm(LBf(interiorv));
    Lmax = max(abs(LBf - LBleg)) / max(abs(LBf));
    Lmaxint = max(abs(LBf(interiorv) - LBleg(interiorv))) / max(abs(LBf(interiorv)));
    fprintf('#P=%d,h=%2.2f Normalized L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', N, emoptions.heath, L2,L2int, Lmax,Lmaxint);
end




%
% euclidean gaussian heat-kernel L-B
%
avgedgelen = edgeStats(mesh.v,mesh.e);
eboptions = [];
eboptions.t = 0.04/4;
Wg = makeEpsBallWeights(P, 'gaussian', 0.3, eboptions);
Ac = vertexArea(mesh,[],'onering') / 3;
WAg = Wg;
for k = 1:n
    WAg(k,:) = Wg(k,:) .* Ac';
end
Dg = sparse(1:n,1:n,sum(WAg,2),n,n);
LBlg = -(Dg-WAg)*f  / (4*pi*eboptions.t^2);

if norm(LBf) == 0
    L2 = norm(LBf-LBlg);
    L2int = norm(LBf(interiorv)-LBlg(interiorv));
    Lmax = max(abs(LBf - LBlg));
    Lmaxint = max(abs(LBf(interiorv) - LBlg(interiorv)));
    fprintf('L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', L2,L2int, Lmax,Lmaxint);
else 
    L2 = norm(LBf - LBlg) / norm(LBf);
    L2int = norm(LBf(interiorv)-LBlg(interiorv)) / norm(LBf(interiorv));
    Lmax = max(abs(LBf - LBlg)) / max(abs(LBf));
    Lmaxint = max(abs(LBf(interiorv) - LBlg(interiorv))) / max(abs(LBf(interiorv)));
    fprintf('#P=%d,h=%2.2f Normalized L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', N, eboptions.t, L2,L2int, Lmax,Lmaxint);
end




%
% heat-kernel L-B, use all vertices
%
h = 0.04/4;
AreaT = faceArea(mesh);
LBlg2 = zeros(size(LBf));
for k = 1:n
    if mod(k,100) == 0
        fprintf('%d / %d finished (%3.2f%%)\n', k, n, (k/n)*100);
    end
    w = mesh.v(k,:);
    W = exp( -vmag2(vadd(mesh.v, -w)) / (4*h) );
    F = vadd(f,-f(k));
    WF = W .* F;
    ksum = 0;
    for fi = 1:numel(mesh.fidx)
        ti = mesh.f(fi,:);
        tsum = sum( WF(ti) );
        ksum = ksum + tsum * AreaT(fi)/3;
    end
    LBlg2(k) = ksum / (4*pi*h^2);
end

if norm(LBf) == 0
    L2 = norm(LBf-LBlg2);
    L2int = norm(LBf(interiorv)-LBlg2(interiorv));
    Lmax = max(abs(LBf - LBlg2));
    Lmaxint = max(abs(LBf(interiorv) - LBlg2(interiorv)));
    fprintf('L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', L2,L2int, Lmax,Lmaxint);
else 
    L2 = norm(LBf - LBlg2) / norm(LBf);
    L2int = norm(LBf(interiorv)-LBlg2(interiorv)) / norm(LBf(interiorv));
    Lmax = max(abs(LBf - LBlg2)) / max(abs(LBf));
    Lmaxint = max(abs(LBf(interiorv) - LBlg2(interiorv))) / max(abs(LBf(interiorv)));
    fprintf('Normalized L2: %f  L2int: %f  Lmax: %f  Lmaxint: %f\n', L2,L2int, Lmax,Lmaxint);
end






%% L-B operator tests on sphere

% points on unit sphere
N = 10;
P = randn(N,3);
P = normalize(P);
Pn = P;
scatter3(P(:,1),P(:,2),P(:,3));


% semi-regular sphere meshes
mesh=readMesh('sphere_semireg_500.obj');
mesh=readMesh('sphere_semireg_1000.obj');
mesh=readMesh('sphere_semireg_1500.obj');
mesh=readMesh('sphere_semireg_2000.obj');
mesh=readMesh('sphere_semireg_4000.obj');

% points on unit sphere
mesh.v = normalize(mesh.v);     
mesh.n = mesh.v;
P = mesh.v;
Pn = mesh.n;
N = size(P,1);
%scatter3(P(:,1),P(:,2),P(:,3));


% L-B functions on spheres


[theta,phi,r] = cart2sphY(P(:,1),P(:,2),P(:,3));

% f = x^2 + y^2
f = P(:,1).^2 + P(:,2).^2;
LBf = 4 * cos(phi).^2 - 2 * sin(phi).^2;










% DEM weights w/o using mesh
emoptions = [];
emoptions.weightmode = 'optimal2';
%emoptions.nbrtype = 'geoball';
%emoptions.nbrsize = 0.2;
%emoptions.nbrsize = 0.025*nbrsize;
emoptions.nbrtype = 'k';
emoptions.nbrsize = 10;
%emoptions.nbrsize = nbrsize;
%emoptions.silent = 1;
emoptions.regtol = 10e-3;
[W,err2,Vu,Vv] = makeExpMapWeights(P, Pn, emoptions);
meanW = mean(nonzeros(W(:)));
    
aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
[A,isboundary] = pointArea(P, [], 'uvVoronoi', aoptions);
Adelsum = pointArea(P, [], 'uvDelArea', aoptions);

WA = W;     EM_A = A;
for k = 1:N
    WA(k,:) = W(k,:) .* A';
    EM_A(k) = (abs(W(k,:))>0)*A + A(k);
end
D = sparse(1:N,1:N,sum(W,2),N,N);
DA = sparse(1:N,1:N,sum(WA,2),N,N);

%LBl = -((DA-WA)*f)/2  ./ ( sum(WA,2) ) ./ EM_A / 0.1;
LBl = -((DA-WA)*f)/2   ./ ( sum(WA,2) ) ./ EM_A / meanW;
%LBl = -((DA-WA)*f)/2  ./ ( pi*realR^2*meanW ) ./ (pi*realR^2) / (meanW*2);
%LBl = -((DA-WA)*f)/2  ./ (EM_A * meanW) ./ EM_A / (meanW * 2);
%LBl = -((DA-WA)*f)/2   * 10000 * 2.4; % ./ (4*pi*0.2^2);
%LBl = -((DA-WA)*f)/2   / (mean(EM_A)*mean(A)); % ./ (4*pi*0.2^2);

%LBl = -(D-W)*f;
%LBl = -(D-W)*f ./ (A) * 2;
%LBl = -(D-W)*f ./ (Adelsum.^2) * 2;
%LBl = -(D-W)*f ./ (EM_A.^2) * 2;

L2 = norm(LBf - LBl) / norm(LBf);
Lmax = max(abs(LBf - LBl)) / max(abs(LBf));
fprintf('[tangent] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);












%
% cotan L-B
%
Wcot = makeOneRingWeights(mesh, 'dcp');
Dcot = sparse(1:N,1:N,sum(Wcot,2),N,N); 
Amix = vertexArea(mesh,[],'mixed');
LBlcot = -(Dcot-Wcot)*f ./ (2*Amix);

L2 = norm(LBf - LBlcot) / norm(LBf);
Lmax = max(abs(LBf - LBlcot)) / max(abs(LBf));
fprintf('[cotan] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);







%
% expmap heat-kernel L-B
%
% DEM weights w/o using mesh
emoptions = [];
emoptions.weightmode = 'gaussian';
emoptions.nbrtype = 'geoball';
emoptions.nbrsize = 0.2;
emoptions.heath = 0.04/4;
[Weg,err2,Vu,Vv] = makeExpMapWeights(P, Pn, emoptions);
Aeg = vertexArea(mesh,[],'onering') / 3;
WAeg = Weg;
for k = 1:N
    WAeg(k,:) = Weg(k,:) .* Aeg';
end
Deg = sparse(1:N,1:N,sum(WAeg,2),N,N);
LBleg = -(Deg-WAeg)*f  / (4*pi*emoptions.heath^2);

L2 = norm(LBf - LBleg) / norm(LBf);
Lmax = max(abs(LBf - LBleg)) / max(abs(LBf));
fprintf('[heatTan] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);



% calculate appropriate epsball size for spherical geoball above
epsballsize = 2 * sin( emoptions.nbrsize / 2 )
heatsize = 2 * sin( emoptions.heath / 2 )



%
% euclidean gaussian heat-kernel L-B
%
avgedgelen = edgeStats(mesh.v,mesh.e);
eboptions = [];
eboptions.t = heatsize;
Wg = makeEpsBallWeights(P, 'gaussian', epsballsize * 0.95, eboptions);
Ac = vertexArea(mesh,[],'onering') / 3;
WAg = Wg;
for k = 1:N
    WAg(k,:) = Wg(k,:) .* Ac';
end
Dg = sparse(1:N,1:N,sum(WAg,2),N,N);
LBlg = -(Dg-WAg)*f  / (4*pi*eboptions.t^2);

L2 = norm(LBf - LBlg) / norm(LBf);
Lmax = max(abs(LBf - LBlg)) / max(abs(LBf));
fprintf('[heat3D] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);



%
% heat-kernel L-B, use all vertices
%
h = 0.04/4;
AreaT = faceArea(mesh);
LBlg2 = zeros(size(LBf));
for k = 1:N
    if mod(k,100) == 0
        fprintf('%d / %d finished (%3.2f%%)\n', k, n, (k/N)*100);
    end
    w = mesh.v(k,:);
    W = exp( -vmag2(vadd(mesh.v, -w)) / (4*h) );
    F = vadd(f,-f(k));
    WF = W .* F;
    ksum = 0;
    for fi = 1:numel(mesh.fidx)
        ti = mesh.f(fi,:);
        tsum = sum( WF(ti) );
        ksum = ksum + tsum * AreaT(fi)/3;
    end
    LBlg2(k) = ksum / (4*pi*h^2);
end

L2 = norm(LBf - LBlg2) / norm(LBf);
Lmax = max(abs(LBf - LBlg2)) / max(abs(LBf));
fprintf('[heatAll] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);





%% compare expmap and euclidean heat-kernels


k_nbrhood = 40;

%
% expmap heat-kernel L-B
%
% DEM weights w/o using mesh
emoptions = [];
emoptions.weightmode = 'gaussian';
emoptions.nbrtype = 'k';
emoptions.nbrsize = k_nbrhood;
emoptions.heath = 0.04/4;
[Weg,err2,Vu,Vv] = makeExpMapWeights(P, Pn, emoptions);
Aeg = vertexArea(mesh,[],'onering') / 3;
%aoptions = [];      aoptions.u = Vu;        aoptions.v = Vv;
%Aeg = pointArea(P, [], 'uvVoronoi', aoptions);
WAeg = Weg;
for k = 1:N
    WAeg(k,:) = Weg(k,:) .* Aeg';
end
Deg = sparse(1:N,1:N,sum(WAeg,2),N,N);
LBleg = -(Deg-WAeg)*f  / (4*pi*emoptions.heath^2);

L2 = norm(LBf - LBleg) / norm(LBf);
Lmax = max(abs(LBf - LBleg)) / max(abs(LBf));
fprintf('[heatTan] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);




%
% euclidean gaussian heat-kernel L-B
%
koptions = [];
koptions.t = heatsize;
Wg = makeKNbrWeights(P, 'gaussian', k_nbrhood, koptions);
Ac = vertexArea(mesh,[],'onering') / 3;
WAg = Wg;
for k = 1:N
    WAg(k,:) = Wg(k,:) .* Ac';
end
Dg = sparse(1:N,1:N,sum(WAg,2),N,N);
LBlg = -(Dg-WAg)*f  / (4*pi*eboptions.t^2);

L2 = norm(LBf - LBlg) / norm(LBf);
Lmax = max(abs(LBf - LBlg)) / max(abs(LBf));
fprintf('[heat3D] Normalized L2: %f  Lmax: %f  \n', L2,Lmax);







%% graphs

% expmap

subplot(1,2,1);
tmpv = [mesh.v(:,1),mesh.v(:,2),abs(LBl-LBf)];
tmpmesh = makeMesh(tmpv, mesh.f, mesh.n);
plotMesh(tmpmesh);

surf(reshape(abs(LBl-LBf), sqrt(N), sqrt(N)));


% cotan

subplot(1,2,2);
tmpv = [mesh.v(:,1),mesh.v(:,2),abs(LBlcot-LBf)];
tmpmesh = makeMesh(tmpv, mesh.f, mesh.n);
plotMesh(tmpmesh);

figure;
surf(reshape(abs(LBlcot-LBf), sqrt(N), sqrt(N)));


