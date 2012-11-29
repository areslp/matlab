%
% [TODO]
%   - try meyer02 mixed-area weights on local 3D meshes from uv-delaunay
%
%
%

%% load a sphere mesh

mesh = readMesh('sphere_graphite_300.obj','n');
mesh = readMesh('sphere_graphite_1000.obj','n');
mesh = readMesh('sphere_graphite_2000.obj','n');
mesh = readMesh('sphere_graphite_4000.obj','n');
mesh = readMesh('sphere_graphite_8000.obj','n');

radius = 1;
mesh.v = radius * normalize(mesh.v);   % project onto sphere
N = size(mesh.v,1);

%% evaluate one of the following laplace-beltrami operators on sphere

% f = x^2 + y^2
vtx_u2 = mesh.v(:,1).^2 + mesh.v(:,2).^2;
[theta,phi,r] = cart2sphM(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
lbf_u2 = (1 + 3*cos(2*theta)) / (radius^2);
figure; plotMesh(mesh,'efb',lbf_u2); title('Analytic');

% f = x + y + z
vtx_u2 = mesh.v(:,1) + mesh.v(:,2) + mesh.v(:,3);
[theta,phi,r] = cart2sphM(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
lbf_u2 = -2*( cos(theta) + (cos(phi)+sin(phi)).*sin(theta) ) / (radius^2);
figure; plotMesh(mesh,'efb',lbf_u2); title('Analytic');


% f = x^2 + (y^2 * z)
vtx_u2 = mesh.v(:,1).^2 + (mesh.v(:,2).^2) .* mesh.v(:,3);
[theta,phi,r] = cart2sphM(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
lbf_u2 = ( 1 + cos(theta) + 3*cos(2*theta) + 3*cos(3*theta) + 6*cos(2*phi).*(-1+2*cos(theta)).*(sin(theta).^2)  ) / (2*radius^2);
figure; plotMesh(mesh,'efb',lbf_u2); title('Analytic');

% f = cos(x)*sin(y)
vtx_u2 = cos(mesh.v(:,1)) .* sin(mesh.v(:,2));
[theta,phi,r] = cart2sphM(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
k1=sin(phi).*sin(theta);
k2=cos(phi).*sin(theta);
lbf_u2 = (1/(radius^2))*(2*k2.*sin(k2).*(cos(k1).*k1 + sin(k1)) - cos(k2).*(2*cos(k1).*k1 + (1 + cos(theta).^2).*sin(k1)) );
figure; plotMesh(mesh,'efb',lbf_u2); title('Analytic'); 


% [RMS] this one oscillates quite a bit on the sphere, to the point where
%   on the 1k mesh at least, heat-kernel weights do not converge w/
%   increasing or decreasing h value. Instead there is a minimum.
%   (expmap weights also diverge as radius increases...basically this
%    function is very non-linear, so larger radius makes things worse...)
% [RMS] actually for this one, on 1k mesh it seems like expmap l-b is more
%   accurate than anything else, w/ radius around [2,3]*edgelen_avg
%     (at least w/ my linear-rescaling metric...)
% f = cos(7*x)*sin(5*y)
vtx_u2 = cos(7*mesh.v(:,1)) .* sin(5*mesh.v(:,2));
[theta,phi,r] = cart2sphM(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
k1=sin(phi).*sin(theta);
k2=cos(phi).*sin(theta);
lbf_u2 = (1/(radius^2))*(14*k2.*sin(7*k2).*(5*cos(5*k1).*k1 + sin(5*k1)) - 0.5*cos(7*k2).*(20*cos(5*k1).*k1 + (37*(3 + cos(2*theta)) - 24*cos(2*phi).*(sin(theta).^2)).*sin(5*k1)) );
figure; plotMesh(mesh,'efb',lbf_u2); title('Analytic');



%% evaluate different L-B approximations

% this only works if 'Voronoi-Laplace' code from
%   http://www.cse.ohio-state.edu/~luoc/Voronoi-Laplace.zip
% is in path 

tangentplane_pts = 15;
ptrtree = BuildGLTree(mesh.v);
[kidc,kdist] = KNNSearch(mesh.v, mesh.v, ptrtree, tangentplane_pts);
DeleteGLTree(ptrtree);
kidc = int32(kidc);
h = (sum(sum(kdist(:, end-10:end-1))) / 10 / N) ^ 2;  
%h = h*2;
%h = h / 5;
AW = VoronoiArea(mesh.v, kidc, sqrt(h));
%AW = ones(size(AW));
G = Gnl_Eigen_Matrix(mesh.v, h, AW);
Wpts = sparse(1 : N, 1 : N, AW);

lbf_E = -G * Wpts * vtx_u2; 
errE_l2 = norm(lbf_E-lbf_u2) / norm(lbf_u2);
errE_linf = norm(lbf_E-lbf_u2,inf) / norm(lbf_u2,inf);
fprintf('error for lbfE:  L2: %8.6f   Linf: %8.6f\n', errE_l2, errE_linf);
figure; plotMesh(mesh,'efb',lbf_E); title(['Euclidean HeatKernel  L2:', num2str(errE_l2), '  Linf:',num2str(errE_linf)]); 


% this one only works if my sphere-geodesics version 'sphereGnl_Eigen_Matrix.m' is in path
AWem = VoronoiArea(mesh.v, kidc, sqrt(h));
Gem = sphereGnl_Eigen_Matrix(mesh.v, h, AWem);
Wptsem = sparse(1 : N, 1 : N, AWem);

lbf_G = -Gem * Wptsem * vtx_u2; 
figure; plotMesh(mesh,'efb',lbf_G); title('Geodesic HeatKernel'); 
errG_l2 = norm(lbf_G-lbf_u2) / norm(lbf_u2);
errG_linf = norm(lbf_G-lbf_u2,inf) / norm(lbf_u2,inf);
fprintf('error for lbfG:  L2: %8.6f   Linf: %8.6f\n', errG_l2, errG_linf);


% cotan w/ voronoi areas
Lcot = - cotanWeights(mesh) / 2;        % cotanWeights function does not div-by-2
Lcot(1:N+1:N*N) = -sum(Lcot,2);  % set W(i,i) = -sum_j W(i,j)
Avor = vertexArea(mesh, [], 'mixed');
for k=1:N
    Lcot(k,:) = Lcot(k,:) / Avor(k);
end

lbf_cot = -Lcot * vtx_u2; 
errCot_l2 = norm(lbf_cot-lbf_u2) / norm(lbf_u2);
errCot_linf = norm(lbf_cot-lbf_u2,inf) / norm(lbf_u2,inf);
fprintf('error for lbfcot:  L2: %8.6f   Linf: %8.6f\n', errCot_l2, errCot_linf);
figure; plotMesh(mesh,'efb',lbf_cot); title(['Cotan w/ Mixed VtxArea:  L2:', num2str(errCot_l2), '  Linf:',num2str(errCot_linf)]); 


% expmap weights
%
%  - with normalization disabled, nbrhoodsize appears to have no
%    significant effect on 'scale' of approx l-b values (ie they
%    always come out w/ roughly the same value)
%       (same appears to go for sphere radius...)
%  - changing regtol changes magnitudes, but otherwise gives same behavior
%    (up to a point...>> regtol means more recons error...)
%
edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
%options.regtol=10e-4;
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;
%options.nbrtype = 'k';
%options.nbrsize = 60;
options.normalize = 0;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
fprintf('recons err is %f\n',sum(err));
Wtan = -Wtan;
Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)

lbf_Wtan = -Wtan * vtx_u2;

errW_l2 = norm(lbf_Wtan-lbf_u2) / norm(lbf_u2);
errW_linf = norm(lbf_Wtan-lbf_u2,inf) / norm(lbf_u2,inf);
fprintf('error for lbfW:  L2: %8.6f   Linf: %8.6f\n', errW_l2, errW_linf);
figure; plotMesh(mesh,'efb',lbf_Wtan); title(['UnNormalized ExpMap:  L2:', num2str(errW_l2), '  Linf:',num2str(errW_linf)]); 



%% area-weighted expmap weights

edgelen_avg = edgeStats(mesh.v, mesh.e);
options = [];
%options.regtol=10e-3;
options.weightmode = 'optimal2';
options.nbrtype = 'geoball';
options.nbrsize = 3.0*edgelen_avg;

% disabling normalization gets us close to correct scale for 300-vtx (but why??)
%   - intrinsic scaling in local weights??
myopts = options;
myopts.normalize = 0;
[Wtan,err,vU,vV] = makeExpMapWeights(mesh, myopts);
ptareaopt = struct('u',vU,'v',vV);

%AW = ones(N,1);
%AW = ones(N,1) /   (mean(abs(lbf_Wtan)) / mean(abs(lbf_u2)));
AW = VoronoiArea(mesh.v, kidc, sqrt(h));
%[AW,isb] = pointArea(mesh.v, [], 'uvVoronoi', ptareaopt);
%[AW,isb] = pointArea(mesh.v, [], 'uvDelRing', ptareaopt);  AW = AW/3;
%[AW,isb] = pointArea(mesh.v, [], 'uvDelRing3', ptareaopt);  AW = AW/3;
%[AW,isb] = pointArea(mesh.v, [], 'uvDelArea', ptareaopt); AW = AW*5;

Lb = -Wtan;
%Lb = Lb * (pi*options.nbrsize^2);
for k = 1:N
    % sum areas of nbr points to get area of ROI disc.... (?)
    % then divide by area of this vertex (??)
    Lb(k,k) = -sum(Lb(k,:) .* AW') / AW(k);
end
Lb = Lb * sparse(1 : N, 1 : N, AW);
% for k = 1:N
%     Lb(k,k) = -sum(Lb(k,:)) / AW(k);
% end

lbf_Wtan = -Lb * vtx_u2;
errW_l2 = norm(lbf_Wtan-lbf_u2) / norm(lbf_u2);
errW_linf = norm(lbf_Wtan-lbf_u2,inf) / norm(lbf_u2,inf);
fprintf('error for lbfW:  L2: %8.6f   Linf: %8.6f\n', errW_l2, errW_linf);

figure; plotMesh(mesh,'efb',lbf_Wtan); title(['AreaWeighted ExpMap:  L2:', num2str(errW_l2), '  Linf:',num2str(errW_linf)]); 



[lbf_Wtan,lbf_u2,lbf_u2./lbf_Wtan]

scale_err = lbf_u2./lbf_Wtan;
[lbf_Wtan,lbf_u2,scale_err]

figure; hold all;
plot(lbf_u2./lbf_Wtan);
plot(AW);
hold off; drawnow;










%%

for k = 2:5

    % expmap weights
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = k*edgelen_avg;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    fprintf('recons err is %f\n',sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)

    lbf_Wtan = -Wtan * vtx_u2;
    WtanS = mean(abs(lbf_Wtan)) / mean(abs(lbf_u2));
    lbf_Wtan = lbf_Wtan/WtanS;
    figure; title('expmap'); plotMesh(mesh,'efb',lbf_Wtan);

    errW_l2 = norm(lbf_Wtan-lbf_u2) / norm(lbf_u2);
    errW_linf = norm(lbf_Wtan-lbf_u2,inf) / norm(lbf_u2,inf);
    fprintf('[%d] error for lbfW:  L2: %8.6f   Linf: %8.6f\n', k, errW_l2, errW_linf);

end




fixed_nbrhood_size = 3 * 0.0822;  % this 3*edgelen_avg for hemisphere_graphite_1000.obj
%meshes = {'hemisphere_graphite_1000.obj', 'hemisphere_graphite_2000.obj', 'hemisphere_graphite_4000.obj', 'hemisphere_graphite_16000.obj'};
meshes = {'hemisphere_graphite_1000.obj', 'hemisphere_graphite_2000.obj', 'hemisphere_graphite_4000.obj'};
for mi = 1:numel(meshes)
    meshname = meshes{mi};
    mesh = readMesh(meshname, 'n');
    mesh.v = normalize(mesh.v);
    
    N = size(mesh.v,1);
    vu2 = mesh.v(:,1).^2 + mesh.v(:,2).^2;
    [theta,phi,r] = cart2sphY(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
    lbfu2 = 1 + 3*cos(2*phi);    
    
    % expmap weights
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
%    options.nbrsize = 3*edgelen_avg;
    options.nbrsize = fixed_nbrhood_size;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    fprintf('recons err is %f\n',sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)

    lbf_Wtan = -Wtan * vu2;
    WtanS = mean(abs(lbf_Wtan)) / mean(abs(lbfu2));
    lbf_Wtan = lbf_Wtan/WtanS;
    figure; title('expmap'); plotMesh(mesh,'efb',lbf_Wtan);

    errW_l2 = norm(lbf_Wtan-lbfu2) / norm(lbfu2);
    errW_linf = norm(lbf_Wtan-lbfu2,inf) / norm(lbfu2,inf);
    fprintf('[%s] error for lbfW:  L2: %8.6f   Linf: %8.6f\n', meshes{mi}, errW_l2, errW_linf);
end


%% testing

ptareaopt = struct('u',vU,'v',vV);
[tanA,isb] = pointArea(mesh.v, [], 'uvVoronoi', ptareaopt);
invD = sparse(1:N,1:N, 1./tanA);

lbft = -Wtan * vtx_u2; 
%lbft = lbft / mean(tanA);
%lbft = lbft / WtanS;
S1 = norm(lbft);  S2 = norm(lbf_u2);
figure; title('expmap'); plotMesh(mesh,'efb',lbft);
errl2 = norm(lbft/S1-lbf_u2/S2) / norm(lbf_u2/S2);
errlinf = norm(lbft/S1-lbf_u2/S2,inf) / norm(lbf_u2/S2,inf);
fprintf('error for lbfW:  L2: %8.6f   Linf: %8.6f\n', errl2, errlinf);


S1 = norm(lbf_E);  S2 = norm(lbf_u2);
errE_l2 = norm(lbf_E/S1-lbf_u2/S2) / norm(lbf_u2/S2);
errE_linf = norm(lbf_E/S1-lbf_u2/S2,inf) / norm(lbf_u2/S2,inf);
fprintf('error for lbfE:  L2: %8.6f   Linf: %8.6f\n', errE_l2, errE_linf);



%%
% rescale all values to linear range (ie color map) and compute error
% there... in this case error seems to decrease as k-nbrhood size
% increases...
exact_range = (lbf_u2-min(lbf_u2))/(max(lbf_u2)-min(lbf_u2));
cot_range = (lbf_cot-min(lbf_cot))/(max(lbf_cot)-min(lbf_cot));
tan_range = (lbf_Wtan-min(lbf_Wtan))/(max(lbf_Wtan)-min(lbf_Wtan));


cot_L2 = errnorm(cot_range,exact_range);
cot_Linf = errnorm(cot_range,exact_range,inf);
[cot_L2,cot_Linf]
tan_L2 = errnorm(tan_range,exact_range);
tan_Linf = errnorm(tan_range,exact_range,inf);
[cot_L2,cot_Linf,tan_L2,tan_Linf]
[cot_L2,cot_Linf]



%% test convergence in color-map for increasing geodesic nbrhood

% This does seem to be converging for increasing nbrhood size.
% But it converges to roughly the same error values for each mesh resolution.

edgelen_avg = edgeStats(mesh.v, mesh.e);

% this is from the 300-vtx sphere...
%edgelen_avg = 0.2072;

do_normalize = 1;       % normalizing seems to have somewhat higher error...
for k = 2:7

    % expmap weights
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = k*edgelen_avg;
    options.normalize = do_normalize;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    fprintf('recons err is %f\n',sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)

    lbf_Wtan = -Wtan * vtx_u2;
    
    exact_range = (lbf_u2-min(lbf_u2))/(max(lbf_u2)-min(lbf_u2));
    tan_range = (lbf_Wtan-min(lbf_Wtan))/(max(lbf_Wtan)-min(lbf_Wtan));
    tan_L2 = errnorm(tan_range,exact_range);
    tan_Linf = errnorm(tan_range,exact_range,inf);    

    figure; plotMesh(mesh,'efb',lbf_Wtan);
    title(['Expmap L-B (colormap errors)  L2:',num2str(tan_L2),'  Linf:',num2str(tan_Linf)]); 
    fprintf('[%d] error for lbfW (r=%8.6f):  L2: %8.6f   Linf: %8.6f\n', k, options.nbrsize, tan_L2, tan_Linf);
end












%% 
% test color-map convergence for different mesh resolutions with fixed
% geodesic nbrhood radius 
%
% 

meshes = {'sphere_graphite_300.obj', 'sphere_graphite_1000.obj', 'sphere_graphite_2000.obj', 'sphere_graphite_4000.obj', 'sphere_graphite_8000.obj'};

fixed_radius = 0.3;
do_normalize = 0;       % normalizing seems to have somewhat higher error...

for mi = 1:numel(meshes)
    mesh = readMesh( meshes{mi}, 'n' );

    radius = 1;
    mesh.v = radius * normalize(mesh.v);   % project onto sphere
    N = size(mesh.v,1);

    vtx_u2 = mesh.v(:,1).^2 + mesh.v(:,2).^2;
    [theta,phi,r] = cart2sphY(mesh.v(:,1), mesh.v(:,2), mesh.v(:,3));
    lbf_u2 = (1 + 3*cos(2*phi)) / (radius^2);

%    figure; plotMesh(mesh,'efb',lbf_u2); title('Analytic');


    % expmap weights
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = fixed_radius;
    options.normalize = do_normalize;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    fprintf('recons err is %f\n',sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)

    lbf_Wtan = -Wtan * vtx_u2;
    
    exact_range = (lbf_u2-min(lbf_u2))/(max(lbf_u2)-min(lbf_u2));
    tan_range = (lbf_Wtan-min(lbf_Wtan))/(max(lbf_Wtan)-min(lbf_Wtan));
    tan_L2 = errnorm(tan_range,exact_range);
    tan_Linf = errnorm(tan_range,exact_range,inf);    

    figure; plotMesh(mesh,'efb',lbf_Wtan);
    title(['Expmap L-B (colormap errors)  L2:',num2str(tan_L2),'  Linf:',num2str(tan_Linf)]); 
    fprintf('[%s] error for lbfW:  L2: %8.6f   Linf: %8.6f\n', meshes{mi}, tan_L2, tan_Linf);
end
    
    


%% finer-grained test for convergence in color-map for increasing geodesic nbrhood

% This does seem to be converging for increasing nbrhood size.
% But it converges to roughly the same error values for each mesh resolution.

edgelen_avg = edgeStats(mesh.v, mesh.e);

% this is from the 300-vtx sphere...
%edgelen_avg = 0.2072;

do_normalize = 0;       % normalizing seems to have somewhat higher error...
for rad = 2*edgelen_avg:edgelen_avg/4:9*edgelen_avg

    % expmap weights
    options = [];
    %options.regtol=10e-5;
    options.weightmode = 'optimal2';
    options.nbrtype = 'geoball';
    options.nbrsize = rad;
    options.normalize = do_normalize;
    [Wtan,err,vU,vV] = makeExpMapWeights(mesh, options);
    fprintf('recons err is %f\n',sum(err));
    Wtan = -Wtan;
    Wtan(1:N+1:N*N) = -sum(Wtan,2);  % set W(i,i) = -sum_j W(i,j)

    lbf_Wtan = -Wtan * vtx_u2;
    
    exact_range = (lbf_u2-min(lbf_u2))/(max(lbf_u2)-min(lbf_u2));
    tan_range = (lbf_Wtan-min(lbf_Wtan))/(max(lbf_Wtan)-min(lbf_Wtan));
    tan_L2 = errnorm(tan_range,exact_range);
    tan_Linf = errnorm(tan_range,exact_range,inf);    

    figure; plotMesh(mesh,'efb',lbf_Wtan);
    title(['Expmap L-B (colormap errors)  L2:',num2str(tan_L2),'  Linf:',num2str(tan_Linf)]); 
    fprintf('[%d] error for lbfW (r=%8.6f):  L2: %8.6f   Linf: %8.6f\n', k, options.nbrsize, tan_L2, tan_Linf);
end
