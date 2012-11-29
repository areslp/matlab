
%% read input meshes

mesh = readMesh('patch.obj');
mesh = readMesh('patch2.obj', 'n');
mesh = readMesh('doghead.obj');
mesh = readMesh('bump_ref.obj');
mesh = readMesh('bump_mc.obj');

mesh = readMesh('doghead_cut.obj');
mesh = readMesh('doghead_cut2.obj');


mesh = readMesh('flatc.obj');
mesh = readMesh('flatc2.obj');
mesh = readMesh('flats.obj');
mesh = readMesh('flat_swissroll.obj');

mesh = readMesh('4bump.obj');

mesh = clipEars(mesh);

plotMesh(mesh,'efb');
% draw with vertices and normals - very slow
%plotMesh(mesh,'vefbn');


%% basic embedding tests

boundaryUV = embedBoundary( mesh, 'circle' );

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.0*edgelen_avg;
weights = makeExpMapWeights(mesh, options);

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 30;
weights = makeExpMapWeights(mesh, options);

mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

tmpmesh = mesh;
tmpmesh.u = embedInterior(tmpmesh, boundaryUV, weights);

subplot(2,1,1);
plotMesh(mesh, 'uefb');
subplot(2,1,2);
plotMesh(tmpmesh, 'uefb');

[ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 1);
fprintf('DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', Dangle, Darea, L2, Linf);    



% contract boundary using weights
mesh2 = mesh;
for i = 1:size(mesh.v,1);
   if mesh.isboundaryv(i)
       nbrs = find(weights(i,:) ~= 0);
       newu = vdot( mesh.u(nbrs,:)', weights(i,nbrs) );
       mesh2.u(i,:) = newu;
   end
end
plotMesh(mesh2, 'uefb');
bi = boundaryUV(:,1);
boundaryUV = [bi, mesh2.u(bi,1), mesh2.u(bi,2)];


% uniform one-ring weights
weights = makeOneRingWeights(mesh, 'uniform');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

% discrete conformal weights
weights = makeOneRingWeights(mesh, 'dcp');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

% discrete authalic weights
weights = makeOneRingWeights(mesh, 'dap');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');

% natural discrete conformal
mesh.u = embedDNCP(mesh);
plotMesh(mesh, 'uefb');

% spectral conformal
mesh.u = embedSCP(mesh);
plotMesh(mesh, 'uefb');



% save parameterized mesh
writeMesh(mesh,'temp.obj');

% save flattened mesh
flatmesh = mesh;
flatmesh.v = [mesh.u(:,1), mesh.u(:,2), zeros(size(mesh.u(:,1)),1)];
writeMesh(flatmesh,'tmpflat.obj');



[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 0, 1);

tmp = mesh;
tmp.v = [mesh.u(:,1), mesh.u(:,2), L2];
plotMesh(tmp);
axis auto;

tmp = mesh;
tmp.v = [mesh.u(:,1), mesh.u(:,2), Cangle];
plotMesh(tmp);
axis auto;

tmp = mesh;
tmp.v = [mesh.u(:,1), mesh.u(:,2), Carea];
plotMesh(tmp);


%% change in distortion as neighbourhood size grows


boundaryUV = embedBoundary( mesh, 'circle' );

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
E = [];
saveU = [];   saveK = [];
for k = 5:5:50
    options.nbrsize = k;
    weights = makeExpMapWeights(mesh, options);
    mesh.u = embedInterior(mesh, boundaryUV, weights);
    [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 1);
    E = [E; Dangle, Darea, L2, Linf ];
    saveU = [saveU, mesh.u];
    saveK = [saveK, k];
    fprintf('k=%d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', k, Dangle, Darea, L2, Linf);    
end


edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
E = [];
saveU = [];   saveK = [];
for r = 1.25:0.25:3.0
    options.nbrsize = r * edgelen_avg;
    weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);
    mesh.u = embedInterior(mesh, boundaryUV, weights);
    [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 1);
    E = [E; Dangle, Darea, L2, Linf ];
    saveU = [saveU, mesh.u];
    saveK = [saveK, r];
    fprintf('r=%5.3f  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', r, Dangle, Darea, L2, Linf);    
end


plot(E(:,1))   % Dangle
plot(E(:,2))   % Darea
plot(E(:,3))   % L2
plot(E(:,4))   % Linf


% plot array of paramterizations
newplot; hold all;
skip = 1;
ucount = size(saveU,2) / 2;
plotsx = ceil(sqrt(ucount / skip));
plotsy = ceil(ucount/skip / plotsx);
for k = 1:skip:ucount
    mesh.u = [ saveU(:,2*k-1), saveU(:,2*k) ];
    subplot(plotsx,plotsy,k/skip)
    plotMesh(mesh, 'uefb');
    title(['param=', num2str(saveK(k))]);
end
hold off; drawnow;
    

% plot sequence of paramterizations  (press spacebar to advance)
skip = 1;
ucount = size(saveU,2) / 2;
plotsx = ceil(sqrt(ucount / skip));
plotsy = ceil(ucount/skip / plotsx);
for k = 1:skip:ucount
    mesh.u = [ saveU(:,2*k-1), saveU(:,2*k) ];
    plotMesh(mesh, 'uefb');
    title(['param=', num2str(saveK(k))]);
    pause;
end




plotMesh(mesh, 'uefb');





%% use LEE embedding

% [RMS TODO] currently manually correcting orientation - should add a
% function for this!!  (cannot avoid w/ LLE)


options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 10;
options.regtol = 1e-5;
weights = makeExpMapWeights(mesh, options);
mesh.u = embedLLE(mesh.v, weights);
plotMesh(mesh, 'uefb');

options.weightmode = 'optimal2';
options.nbrtype = 'k';
E = [];
for k = 5:3:30
    options.nbrsize = k;
    options.regtol = 1e-1;
    weights = makeExpMapWeights(mesh, options);
    mesh.u = embedLLE(mesh.v, weights);
    [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 1);
    if Darea < 0
        mesh.u(:,1) = -mesh.u(:,1);     % correct orientation
        [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 1);
    end
    E = [E; Dangle, Darea, L2, Linf ];
    fprintf('k=%d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', k, Dangle, Darea, L2, Linf); 
    
    newplot;
    plotMesh(mesh, 'uefb');
    drawnow;
end




%% curvature-scaled nbrhoods test

boundaryUV = embedBoundary( mesh, 'circle' );

% find mean-curvature weights
vH = meanCurv(mesh);
vH = abs(vH);
minH = min(vH);   maxH = max(vH);
vH = vadd(vH,-minH);
vH = vH / (maxH-minH);
vC = [vH,vH,vH];
plotMesh(mesh,'efb',vC)

vW = ones(size(vH)) - vH;

NSize = 20;

% compute w/ 
edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 10 + NSize;
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);
clear('options');

mesh1 = mesh;
mesh1.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh1, 'uefb');
[L2,Linf,TareaV,TareaU] = triStretch(mesh1.v, mesh1.u, mesh1.f);
L2 = (L2 .* L2 .* TareaV) / sum(sum(TareaV));
t1 = sort(L2);


edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.vnbrsize = ceil( 10 + ones(size(vH)) * NSize .* vW );
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);
clear('options');

mesh2 = mesh;
mesh2.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh2, 'uefb');
[L2,Linf,TareaV,TareaU] = triStretch(mesh2.v, mesh2.u, mesh2.f);
L2 = (L2 .* L2 .* TareaV) / sum(sum(TareaV));
t1 = sort(L2);



edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.vnbrsize = ceil( 10 + ones(size(vH)) * NSize .* vH );
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);
clear('options');

mesh3 = mesh;
mesh3.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh3, 'uefb');
[L2,Linf,TareaV,TareaU] = triStretch(mesh3.v, mesh3.u, mesh3.f);
L2 = (L2 .* L2 .* TareaV) / sum(sum(TareaV));
t1 = sort(L2);


n = size(t1,1);
%nskip = ceil(n/2);
nskip = 50;
newplot; hold all;
plot(t1(1:n-nskip,:));  plot(t2(1:n-nskip,:));  plot(t3(1:n-nskip,:));
hold off; drawnow;


%% plot graphs of sorted per-vertex weights
%  (kind of interesting patterns)


interiorv=mesh.vidx(mesh.isboundaryv==0);
boundaryv=mesh.vidx(mesh.isboundaryv==1);
filterv = interiorv;
filterv = boundaryv;
filterv = mesh.vidx;
filterv = [1];

% plot graphs for varying 
reglow=3;
reghi=3;
nplots=ceil(sqrt(reghi-reglow+1));
figure;
hold all;
for k = reglow:reghi

    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = 3.0*edgelen_avg;
    options.regtol = 10^-k;
    [weights,err2,vu,vv] = makeExpMapWeights(mesh.v, mesh.n, options);
    fprintf('approximation error is %f\n',err2);
    geodist = sqrt( vu.^2 + vv.^2 );
    
    % normalize weights
    %geodist = vmul(geodist,1./sum(geodist,2));

    % filter to only interior verts
    weightsf=weights(filterv,:);
    geodistf=geodist(filterv,:);
    
    idx=find(weightsf);
    W = full(weightsf(idx));
    G = full(geodistf(idx));
    subplot(nplots,nplots,k-reglow+1);
    meang=mean(G);
    line([meang,meang],[min(W),max(W)]);
    meanmaxg=mean(max(geodist'));
    line([meanmaxg,meanmaxg],[min(W),max(W)]);
    line([min(G),max(G)],[0,0]);
    scatter(G,W);
end
hold off;
drawnow;




figure;
boundaryUV = embedBoundary( mesh, 'circle' );
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');



% [TODO] use this to experiment with regularization constant

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 15;
[weights,err2] = makeExpMapWeights(mesh.v, mesh.n, options, mesh);

hold all;
for i = 1:size(mesh.v,1)
     if mesh.isboundaryv(i)
         continue;
     end
    nzw = sort( weights(i, weights(i,:) > 0 ) );
    plot(nzw)
end
hold off;
drawnow;


    
%% search for optimal regularization tolerance

% results are very insensitive to value... anything in range
% 1e-1 to 1e-8 gives ok results. < 1e-4, nothing much changes...

boundaryUV = embedBoundary( mesh, 'circle' );

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 10;
E = [];
for k = 0.5:0.25:5
    options.regtol = 10^-k;
    [weights,err2] = makeExpMapWeights(mesh.v, mesh.n, options, mesh);
    mesh.u = embedInterior(mesh, boundaryUV, weights);
    [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 1);
    E = [E; k, Dangle, Darea, L2, Linf, err2 ];
    fprintf('k=%4.2f  DAngle: %f  Darea: %f  L2: %f  Linf: %f  Rerr: %f\n', k, Dangle, Darea, L2, Linf, err2);    
end

plot(E(:,1), E(:,2))
plot(E(:,1), E(:,3))
plot(E(:,1), E(:,4))
plot(E(:,1), E(:,5))
plot(E(:,1), E(:,6))

newplot; hold all;
plot(E(:,1), E(:,2) / E(1,2))
plot(E(:,1), E(:,3) / E(1,3))
plot(E(:,1), E(:,4) / E(1,4))
plot(E(:,1), E(:,5) / E(1,5))
hold off; drawnow;


newplot; hold all;
plot(E(:,1), E(:,2) / min(E(:,2)))
plot(E(:,1), E(:,3) / min(E(:,3)))
plot(E(:,1), E(:,4) / min(E(:,4)))
plot(E(:,1), E(:,5) / min(E(:,5)))
hold off; drawnow;




%% linear position constraints

% I thought constraints might work better with large
% nbrhoods, but it doesn't seem like it - the same 
% foldovers occur. I guess maybe this makes sense - larger
% nbrhoods pull the parameterization away from the 
% constrained boundary, so why not away from interior
% fixed points as well? 
%
% (what about anisotropic nbrhood sizes?)

% uniform one-ring weights
weights = makeOneRingWeights(mesh, 'uniform');

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 2.5*edgelen_avg;
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 10;
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);


mesh.u = embedInterior(mesh, boundaryUV, weights);
cpt = [0.25,0];    % point closest to this one gets constrained
dists = vmag2(vadd(mesh.u,-cpt));
[minv,mini]=min(dists);
constraints = [mini,0.35,0];

mesh.u = embedInterior(mesh, boundaryUV, weights);
subplot(1,2,1);
plotMesh(mesh, 'uefb');
mesh.u = embedInterior(mesh, boundaryUV, weights, constraints);
subplot(1,2,2);
plotMesh(mesh, 'uefb');







%% iterative boundary update

% THIS WORKS! W00T!!
%   - only seems to converge for flat meshes...otherwise we see 
%     a u-shape for most of the distortion metrics...
%   - cannot 'expand' out compressed areas of mesh (basically can
%       only move the boundary inwards)
%   - contractBoundary() has hacks to try to deal with this, but
%       it doesn't really work...
%   - what about 'edge-tweaking' from virtual boundary paper?
%

tmpmesh = mesh;
tmpmesh.v = [mesh.v(:,1), mesh.v(:,3), mesh.v(:,2)];
plotMesh(tmpmesh,'efbl');

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 15;
weights = makeExpMapWeights(mesh.v, mesh.n, options);

% geodist seems to work better than k-param for this...
edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.5*edgelen_avg;
weights = makeExpMapWeights(mesh.v, mesh.n, options);

% doesn't work with this...
weights = makeOneRingWeights(mesh, 'uniform');
weights = makeOneRingWeights(mesh, 'invdist2');

plotMesh(mesh, 'efb');

boundaryUV = embedBoundary( mesh, 'circle' );
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');


% iterate 1) b = W b  2) solve Mu = [0 b]
compute_error = 0;
boundaryUV = embedBoundary( mesh, 'circle' );
E = [];
for ITER = 1:100
    
    mesh.u = embedInterior(mesh, boundaryUV, weights);

    if ITER == 1
        [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 0);
        E = [E; Dangle, Darea, L2, Linf ];    
        fprintf('iter %d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', ITER, Dangle, Darea, L2, Linf);    
    end      
    
    % reconstruct boundary from weights
    %  (current boundary is used...)
    for ITER2 = 1:1
        updateu = mesh.u;
        for i = 1:size(mesh.v,1);
           if mesh.isboundaryv(i)
               nbrs = find(weights(i,:) ~= 0);
               newu = vdot( mesh.u(nbrs,:)', weights(i,nbrs) );
               updateu(i,:) = newu;
           end
        end
        mesh.u = updateu;
    end
    bi = boundaryUV(:,1);
    boundaryUV = [bi, mesh.u(bi,1), mesh.u(bi,2)];

    newplot;
    plotMesh(mesh, 'uefb');
    drawnow;    
    
    if compute_error     
        [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 0);
        E = [E; Dangle, Darea, L2, Linf ];    
        fprintf('iter %d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', ITER, Dangle, Darea, L2, Linf);    
    end
end



% iterate u = W u
compute_error = 0;
boundaryUV = embedBoundary( mesh, 'circle' );
mesh.u = embedInterior(mesh, boundaryUV, weights);
E = [];
for ITER = 1:100
    if ITER == 1
        [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 0);
        E = [E; Dangle, Darea, L2, Linf ];    
        fprintf('iter %d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', ITER, Dangle, Darea, L2, Linf);    
    end       
    
    mesh.u = [weights*mesh.u(:,1), weights*mesh.u(:,2)];
%    c = norm( [weights*mesh.u(:,1), weights*mesh.u(:,2)] );    
%    mesh.u = mesh.u / c;    % scale by approx eigenvalue magnitude (as in iterative power method)

    mesh.u = vadd(mesh.u,-mean(mesh.u));
    [V,D] = eigs(cov(mesh.u));    
%     dscale = D(1,1) / D(2,2);
%     e1 = normalize(V(:,1)');   e2 = normalize(V(:,2)');
%     uproj = vdot(mesh.u,e2);
%     ushift = dscale*[uproj * e2(1), uproj * e2(2)];
%     mesh.u = mesh.u + ushift;
        
    newplot;
    hold all;
    plotMesh(mesh, 'uefb');
    drawline([0,0],D(1,1)*V(:,1)')
    drawline([0,0],D(2,2)*V(:,2)')
    hold off;
    %fprintf('c is %f,  D1 %f, D2 %f\n', c, D(1,1), D(2,2));
    drawnow;  

    if compute_error 
        [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 0);
        E = [E; Dangle, Darea, L2, Linf ];    
        fprintf('iter %d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', ITER, Dangle, Darea, L2, Linf);    
    end
end


% play around with eigenvectors of weight matrix...

[V,D] = eigs(weights);

for k = 1:size(V,2)
    fprintf('%d:  %f   %f\n', k, dot(V(:,k), mesh.u(:,1)), dot(V(:,k), mesh.u(:,2)));
end

tmp = mesh;
tmp.u = [V(:,4), V(:,5)];
plotMesh(tmp, 'uefb');

N = size(weights,1);
W = [weights, zeros(size(weights)); zeros(size(weights)), weights ];
[V,D] = eigs(W);
tmp = mesh;
useV = 2;
tmp.u = [V(1:N,useV), V(N+1:2*N,useV)];
plotMesh(tmp, 'uefb');



% iterate u = (1-alpha)*W u + (alpha)*u
alpha = 0.1;
beta = 0.0;
boundaryUV = embedBoundary( mesh, 'circle' );
mesh.u = embedInterior(mesh, boundaryUV, weights);
bi = mesh.loops;
binext = circshift(bi,1);
biprev = circshift(bi,-1);
vup = normalize(mesh.u);
for ITER = 1:250
    newu = [weights*mesh.u(:,1), weights*mesh.u(:,2)];
    mesh.u = (1-alpha)*mesh.u + (alpha)*newu;
    ucentroid = 0.5 * ( mesh.u(binext,:) + mesh.u(biprev,:) );
    ulaplacian = mesh.u(bi,:) - ucentroid;
    ulaplacian = normalize(ulaplacian);
    mesh.u(bi,:) = mesh.u(bi,:) + beta*ulaplacian;
    newplot;
    plotMesh(mesh, 'uefb');
    drawnow;  
end



for i = 1:4
    subplot(2,2,i)
    plot(E(:,i))
end


for i = 1:4
    subplot(2,2,i)
    plot(E(1:25,i))
end





% try to correct for nbrhood contraction by
% smoothing and expanding boundary (hack!)
boundaryUV = embedBoundary( mesh, 'circle' );
E = [];
for ITER = 1:100

    mesh.u = embedInterior(mesh, boundaryUV, weights);
    updateu = contractBoundary( mesh, boundaryUV, weights );
    
    mesh.u = updateu;
    bi = boundaryUV(:,1);
    boundaryUV = [bi, mesh.u(bi,1), mesh.u(bi,2)];
  
    newplot;
    plotMesh(mesh, 'uefb');
    drawnow;    
    
%    [ Dangle, Darea, L2, Linf ] = meshDistortion(mesh, 0);
%    E = [E; Dangle, Darea, L2, Linf ];    
%    fprintf('iter %d  DAngle: %f  Darea: %f  L2: %f  Linf: %f\n', ITER, Dangle, Darea, L2, Linf);    
end


% save parameterized mesh
writeMesh(mesh,'temp.obj');


%% direct boundary embedding

vbi = find(mesh.isboundaryv == 1);
vb = lle(mesh.v(vbi,:)', 4, 2)';
boundaryUV = [vbi, vb(:,1), vb(:,2)];



boundaryUV = embedBoundary( mesh, 'circle' );

weights = makeOneRingWeights(mesh, 'uniform');
mesh.u = embedInterior(mesh, boundaryUV, weights);
plotMesh(mesh, 'uefb');





%% estimate mean curvature normals  (Hn)

% - increasing nbrhood size results in smoother H (nice!)
% - H needs to be scaled somehow...(relative values look good in colormap)
% - n randomly points in very wrong direction
%     - maybe this is telling us about some sort of failure in weights??

mesh = readMesh('patch.obj');
mesh = readMesh('doghead.obj');
mesh.v = [mesh.v(:,1), mesh.v(:,3), mesh.v(:,2)];
mesh.n = [mesh.n(:,1), mesh.n(:,3), mesh.n(:,2)];

edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'k';
options.nbrsize = 15;
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);


edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
weights = makeExpMapWeights(mesh.v, mesh.n, options, mesh);


% existing methods
vH = meanCurv(mesh, mesh.vidx, 'normal');
plotMesh(mesh,'fb', abs(vH))
plotMesh(mesh,'fb', vH)

[vH,Hn] = meanCurv(mesh, mesh.vidx, 'dmc');
vH(mesh.isboundaryv~=0) = 0;
vE = 1 - abs(vdot(mesh.n,normalize(Hn)));
plotMesh(mesh,'fb', abs(vH))


% compute Hn using nbrhood weights, then abs(H) = |Hn|
colormap default;
vH = zeros(size(mesh.v,1),1);
vE = zeros(size(mesh.v,1),1);
Hnormals = mesh.n;
for i = 1:size(mesh.v,1);
   nbrs = find(weights(i,:) ~= 0);
   newv = vdot( mesh.v(nbrs,:)', weights(i,nbrs) )';
   Hn = mesh.v(i,:) - newv;
   vH(i) = vmag(Hn);
   newn = normalize(Hn);
   vE(i) = 1 - abs(vdot(newn, mesh.n(i,:)));
   Hnormals(i,:) = newn;
end
plotMesh(mesh,'fb', abs(vH))        % abs(H)


% compute average nbrhood size
nw = mesh.vidx;
for i = 1:size(mesh.v,1)
    nw(i) = nnz(weights(i,:));
end
mean(nw)
    
% plot normal error
colormap summer;
colormap bone;
plotMesh(mesh,'fb', vE)             

tmpmesh = mesh;
tmpmesh.n = normalize(Hn);
plotMesh(tmpmesh,'fbn',vE);

plotMesh(mesh,'fbl');
plotMesh(mesh,'efbn');


tmpmesh = mesh;
tmpmesh.n = normalize(Hn);
plotMesh(tmpmesh,'fbl');


plotMesh(mesh,'fi');





%% tmp

% points on plane
N = 500;
uv = genpoints(N,2,'stratified_unit_square');
N = size(uv,1);
TRI = delaunay(uv(:,1), uv(:,2));
mesh = makeMesh( [uv,zeros(N,1)], TRI, repmat([0,0,1],N,1) );
figure;
plotMesh(mesh);
P = mesh.v;
Pn = mesh.n;
n = size(P,1);
figure;
scatter3(P(:,1),P(:,2),P(:,3));

% L-B functions on planes

f = 2*P(:,1);           %  f = 2x
LBf = zeros(size(f));   %  f_xx + f_yy = 0

f = P(:,1).^2;          %  f = x^2
LBf = 2*ones(size(f));  %  f_xx + f_yy = 2


% points on sphere
N = 500;
P = randn(N,3);
P = normalize(P);
n = size(P,1);
scatter3(P(:,1),P(:,2),P(:,3));

% points on sphere
mesh=readMesh('sphere.obj');
P = mesh.v;
n = size(P,1);
scatter3(P(:,1),P(:,2),P(:,3));


% L-B functions on spheres

f = 0;  % ...



% DEM weights w/o using mesh
emoptions = [];
emoptions.weightmode = 'optimal2';
%emoptions.nbrtype = 'geoball';
%emoptions.nbrsize = 0.3;
emoptions.nbrtype = 'k';
emoptions.nbrsize = 10;
[W,err2,Vu,Vv] = makeExpMapWeights(P, Pn, emoptions);
aoptions = [];
aoptions.u = Vu;
aoptions.v = Vv;
[A,isboundary] = pointArea(P, [], 'uv', aoptions);
WA = W;
EM_A = A;
for k = 1:n
    WA(k,:) = W(k,:) .* A';
    EM_A(k) = (W(1,:)>0)*A + A(k);
end
%D = sparse(1:n,1:n,sum(WA,2),n,n);
%LBl = ((D-WA)*f) ./ (EM_A);
D = sparse(1:n,1:n,sum(W,2),n,n);
LBl = (D-W)*f;

L2 = norm(LBf - LBl) / norm(LBf);
Lmax = max(abs(LBf - LBl)) / max(abs(LBf));
fprintf('L2: %f   Lmax: %f\n', L2,Lmax);





% cotan weights
W = makeOneRingWeights(mesh, 'dcp');
D = sparse(1:n,1:n,sum(W,2),n,n); 
A = vertexArea(mesh,[],'mixed');
LBl = (D-W)*f;

L2 = norm(LBf - LBl) / norm(LBf);
Lmax = max(abs(LBf - LBl)) / max(abs(LBf));
fprintf('L2: %f   Lmax: %f\n', L2,Lmax);


% DEM weights
edgelen_avg = edgeStats(mesh.v, mesh.e);
emoptions = [];
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
[W,err2,Vu,Vv] = makeExpMapWeights(P, Pn, options);
aoptions = [];
aoptions.u = Vu;
aoptions.v = Vv;
[A,isboundary] = pointArea(mesh.v, [], 'uv', aoptions);
WA = W;
for k = 1:n
    WA(k,:) = W(k,:) .* A';
end
D = sparse(1:n,1:n,sum(WA,2),n,n);
LBl = (D-W)*f;

L2 = norm(LBf - LBl) / norm(LBf);
Lmax = max(abs(LBf - LBl)) / max(abs(LBf));
fprintf('L2: %f   Lmax: %f\n', L2,Lmax);


% heat-kernel / gaussian weights
eboptions = [];
eboptions.t = 0.1/4;
W = makeEpsBallWeights(P, 'gaussian', 10*eboptions.t, eboptions);
%A = vertexArea(mesh,[],'onering') / 3;
for k = 1:n
    WA(k,:) = W(k,:) .* A';
end
D = sparse(1:n,1:n,sum(WA,2),n,n);
LBl = (1/(4*pi*eboptions.t^2)) * (D-WA)*f;

L2 = norm(LBf - LBl) / norm(LBf);
Lmax = max(abs(LBf - LBl)) / max(abs(LBf));
fprintf('L2: %f   Lmax: %f\n', L2,Lmax);



