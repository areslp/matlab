function [ mesh ] = fair_LBH( init_mesh, nMaxIters, convergenceMult, show_progress )
%FAIR_LBH generate mesh that is fair wrt L-B of Mean Curvature H
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   This is an implementation of Schneider and Kobbelt 01

if ~exist('show_progress','var')
    show_progress = 0;  % plot mesh at each iteration (slow but fun to watch)
end

if ~exist('nMaxIters','var')
    nMaxIters = 100;       
end

if ~exist('convergenceMult','var')
    convergenceMult = 0.0001;
end

%t_s = 0.9;          % scaling factor for vertex updates - paper says 0.9
t_s = 0.25;          % scaling factor for vertex updates - paper says 0.9

% laplacian defines inner fairness condition
laplacian_mode = 'uniform';
%laplacian_mode = 'cotan';

% convergence when average vtx movement < stop_criteria
% (this value may be a bit conservative...)
meshradius = vmag(init_mesh.bounds(1,:)-init_mesh.bounds(2,:));
stop_criteria = meshradius * convergenceMult;


mesh = init_mesh;
n = size(mesh.v,1);
Bi = mesh.vidx( mesh.isboundaryv == 1 );  % boundary verts 
Ii = mesh.vidx( mesh.isboundaryv == 0 );  % interior verts 

converged = 0;
delta = 0;
for iter = 1:nMaxIters
    startv = mesh.v;
    
    fprintf('iter %d  (last delta: %f   stop_critera: %f\n', iter, delta, stop_criteria);

    % compute current H's using normal curvature method
    tic; fprintf('   ...computing mean curvatures');
    H = meanCurv( mesh, [], 'normal');
    fprintf('\t\t\t\t\t( %f seconds )\n', toc);

    % construct Dirichlet system for curvatures using cotan weights (Equations 6 & 7)
    tic; fprintf('   ...finding new interior H values');
    W = cotanWeights(mesh);
    S = spdiags(sum(W,2),0,n,n) - W;
    S(Bi,:) = 0;
    S = S + spdiags(mesh.isboundaryv,0,n,n);
    b = H;
    b(Ii) = 0;

    % smoothly interpolate boundary H's over interior
    newH = S \ b;
    fprintf('\t\t\t\t( %f seconds )\n', toc);

    % now update each vertex
    tic;fprintf('   ...Updating vertices for new H values');
%    for i = Ii
    for vi = 1:numel(Ii);
        i = Ii(vi);

        [Hi,data] = meanCurv(mesh, i, 'normal');
        qi = mesh.v(i,:);
        ni = -mesh.n(i,:);  % reverse normal (math assumes inward normals)
        Qj = data(:,1:3);
        Tx = data(:,4);   Ty = data(:,5);
        A = [Tx.^2, Tx.*Ty, Ty.^2];
        M = inv(A'*A)*A';
        denom = vmag2( Qj );
        Cj = 2 * vdot( Qj, ni ) ./ denom;  
        Dj = -2 ./ denom;

        Mc = (M(1,:) + M(3,:)) * Cj;
        Md = (M(1,:) + M(3,:)) * Dj;
        t = ( 2*newH(i) - Mc) / Md;

        % compute updated qi (inner fairness condition)
        Li = laplacian(mesh,i, laplacian_mode);
        qk = qi + Li - vdot(Li,ni)*ni;
 
        mesh.v(i,:) = qk + t_s*t*ni;
%        mesh.v(i,:) = qi + t_s*t*ni;   % update w/o inner fairness
    end
    fprintf('\t\t( %f seconds )\n', toc);

    % estimate new interior normals
    tic;fprintf('   ...Estimating new normals');
    mesh.n(Ii,:) = estimateNormal(mesh,Ii,'faceavg');
    fprintf('\t\t\t\t\t( %f seconds )\n', toc);
    
    if show_progress
        plotMesh(mesh,'efbn');
        drawnow;
    end
    
    delta = mean(vmag(mesh.v-startv));
    if ( delta < stop_criteria )
        fprintf('Converged in %d steps!\n', iter);
        converged = 1;
        break;
    end    
     
end

if ~converged
    fprintf('Reached max iterations before convergence...\n');
end
    
end
