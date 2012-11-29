

mesh=readMesh('sphere.obj');
N = numel(mesh.vidx);


interiorv=mesh.vidx(mesh.isboundaryv==0);
boundaryv=mesh.vidx(mesh.isboundaryv==1);
filterv = interiorv;
filterv = boundaryv;
filterv = mesh.vidx;
filterv = [1];



%% plot graphs of sorted per-vertex weights
%  (kind of interesting patterns)




% plot graphs for varying 
reglow=3;
reghi=3;
nplots=ceil(sqrt(reghi-reglow+1));
figure;
hold all;
for k = reglow:reghi

    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options.weightmode = 'optimal2';
%    options.nbrtype = 'geoball';
%    options.nbrsize = 3.0*edgelen_avg;
    options.nbrtype = 'k';
    options.nbrsize = 8;
    options.regtol = 10^-k;
    [weights,err2,vu,vv] = makeExpMapWeights(mesh.v, mesh.n, options);
    fprintf('approximation error is %f\n',mean(err2));
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




figure;
hold on;
geodist = sqrt( Vu.^2 + Vv.^2 );
weightsf=W(filterv,:);
geodistf=geodist(filterv,:);
idx=find(weightsf);
Wf = full(weightsf(idx));
Gf = full(geodistf(idx));
meang=mean(Gf);
line([meang,meang],[min(Wf),max(Wf)]);
meanmaxg=mean(max(geodist'));
line([min(Gf),max(Gf)],[0,0]);
scatter(Gf,Wf);
hold off;
drawnow;






%% distance and angle weight-distributions


edgelen_avg = edgeStats(mesh.v, mesh.e);
options.weightmode = 'optimal2';
%options.nbrtype = 'geoball';
%options.nbrsize = 3.0*edgelen_avg;
options.nbrtype = 'k';
options.nbrsize = 8;
options.regtol = 10^-3;
[Wtan,err2,vu,vv] = makeExpMapWeights(mesh.v, mesh.n, options);
fprintf('approximation error is %f\n',mean(err2));
geodist = sqrt( vu.^2 + vv.^2 );
    

%
% plot weight by point-angle   (each nbr gets half of angle to prev & next
%                               nbrs, when sorted by angle in tangent space)
%
allangles = [];
allweightsA = [];
for i = 1:N
    nbrs = find(Wtan(i,:) > 0);
    K = numel(nbrs);
    u = full(vu(i,nbrs));
    v = full(vv(i,nbrs));
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
    
    weights = full(Wtan(i,nbrs));   
    allangles = [allangles; nbrangles];
    allweightsA = [allweightsA; weights'];
end


figure;
hold all;
line([min(allangles),max(allangles)],[0,0]);
scatter(allangles,allweightsA);
hold off;
drawnow;




%% same for cotangent weights


Wcot = makeOneRingWeights(mesh, 'dcp', 1);


%
% plot weight by distance
%
alldists = [];
allweightsD = [];
for k = 1:N
    nbrs = find(mesh.e(k,:) > 0);
    dists = vmag(vadd(mesh.v(nbrs,:), -mesh.v(k,:)));
    weights = nonzeros(Wcot(k,:));
    alldists = [alldists; dists];
    allweightsD = [allweightsD; weights];
end


figure;
hold all;
line([min(alldists),max(alldists)],[0,0]);
scatter(alldists,allweightsD);
hold off;
drawnow;


%
% plot weight by vertex-angle   (each nbr gets half of angle of each
%                                triangle around vtx)
%

allangles = [];
allweightsA = [];
angleMat = mesh.e;  
angleMat(:) = 0;
for i = 1:N
    vTris = oneringf(mesh, i);
    for ti = 1:numel(vTris)
        [j,k] = tripick(mesh.f(vTris(ti),:), i);
        e1 = normalize( mesh.v(j,:) - mesh.v(i,:) );
        e2 = normalize( mesh.v(k,:) - mesh.v(i,:) );
        t = vangle(e1,e2);
        angleMat(i,j) = angleMat(i,j) + t/2;
        angleMat(i,k) = angleMat(i,k) + t/2;
    end
    nbrs = find(mesh.e(i,:) > 0);
    angles = full(angleMat(i,nbrs));
    weights = full(Wcot(i,nbrs));    
    allangles = [allangles; angles'];
    allweightsA = [allweightsA; weights'];
end


figure;
hold all;
line([min(allangles),max(allangles)],[0,0]);
scatter(allangles,allweightsA);
hold off;
drawnow;


