function [ deformed_mesh ] = meshDeform_Lipman04( mesh, constraints, rotation_iters, smooth_radius )
%MESHDEFORM_LIPMAN04 Laplacian mesh deformation using estimated rotations
%   implementation of Lipman et al SMI04 paper 
%     "Differential Coordinates for Interactive Mesh Editing"
%
%   constraints:      rows of [vtx_i, x, y, z, weight_i]
%   rotation_iters:   number of rotation estimation iterations   (0 gives pure laplacian deformation)
%   smooth_radius:    radius for normal smoothing. no smoothing if 0. estimated if undefined.
%      [TODO] actually support this parameter as described in paper...

deformed_mesh = mesh;

if ~ exist('rotation_iters', 'var')
    rotation_iters = 1;
end
if ~ exist('smooth_radius', 'var')
    smooth_radius = -1;
end



%% part 1 - solve initial Laplacian deformation, without rotation vectors

% set up Laplacian matrix
M = mesh.e(mesh.vidx,:);
nVerts = numel(mesh.vidx);
nbr_i = zeros(nVerts,1);         % pick most-orthogonal nbr for each vertex
for qi = 1:nVerts
    [ot,ov] = onering(mesh, qi);
    nNbrs = numel(ov);

    % need to find most-orthogonal outgoing edge
    n = mesh.n(qi,:);  v = mesh.v(qi,:);
    nbrdots = zeros(nNbrs,1);
    
    for j = 1:numel(ov)
        qj = ov(j);
        
        % uniform weights
        M(qi,qj) = -1/nNbrs;   
        
        % compute dot product w/ normalized vector to nbr
        vj = vadd(mesh.v(qj,:), -v);  
        vj = normalize(vj);
        nbrdots(j) = abs(vdot(vj,n));
    end
    M(qi,qi) = 1;
    
    % pick out most-orthogonal outgoing edge
    [mindot,nbri] = min(nbrdots);
    nbr_i(qi) = ov(nbri);
end

% compute differential coordinates
D = [ M*mesh.v(:,1), M*mesh.v(:,2), M*mesh.v(:,3) ];


% add soft position constraints
nCons = size(constraints,1);
for i = 1:nCons
   ci = constraints(i,1);
   wi = constraints(i,5);
   nrows = size(M,1);
   M(nrows+1, ci) = wi;
   D(nrows+1, :) = wi * constraints(i,2:4);
end

% least-squares solve for deformed roi
deformed_mesh.v(:,1) = M \ D(:,1);
deformed_mesh.v(:,2) = M \ D(:,2);
deformed_mesh.v(:,3) = M \ D(:,3);

deformed_mesh.n = estimateNormal(deformed_mesh);

if rotation_iters == 0
    return;
end


%% Part2 - estimate rotation for each laplacian vector

orig_D = D;
if  smooth_radius ~= 0
    mesh.n = smooth_normals(mesh, smooth_radius);
end

for iter = 1:rotation_iters

    % smooth normals on deformed mesh
    if  smooth_radius ~= 0 
        deformed_mesh.n = smooth_normals(deformed_mesh, smooth_radius);
    end
    
    for qi = 1:nVerts
        qj = nbr_i(qi);

        % find original frame
        pi = mesh.v(qi,:);
        ni = mesh.n(qi,:);
        pj = mesh.v(qj,:);
        uij = vadd(pj,-pi);
        uij = vadd(uij, -vdot(uij,ni)*ni);
        uij = normalize(uij);
        e2 = ncross(ni, uij);

        % transfer differential coords to frame coords
        d = orig_D(qi,:);
        alpha = vdot(ni,d);
        beta = vdot(uij,d);
        gamma = vdot(e2,d);

        % construct new frame at deformed vertex, with estimated normal
        pi = deformed_mesh.v(qi,:);
        ni = deformed_mesh.n(qi,:);
        pj = deformed_mesh.v(qj,:);
        uij = vadd(pj,-pi);
        uij = vadd(uij, -vdot(uij,ni)*ni);
        uij = normalize(uij);
        e2 = ncross(ni, uij);   

        new_d = alpha*ni + beta*uij + gamma*e2;

        D(qi,:) = new_d;
    end


    % least-squares solve for deformed roi
    deformed_mesh.v(:,1) = M \ D(:,1);
    deformed_mesh.v(:,2) = M \ D(:,2);
    deformed_mesh.v(:,3) = M \ D(:,3);
    deformed_mesh.n = estimateNormal(deformed_mesh);
end



end





% average normals at each vertex, within estimated geodesic radius
% [TODO] actually just averages one-ring normals for now...need to
%    integrate dijkstra nbrhood-finding
function [ N ] = smooth_normals( mesh, smooth_radius )
    N = mesh.n;
   
    % [TODO] support this parameter...
    smooth_radius = -1;
    
    nVerts = numel(mesh.vidx);
    for qi = 1:nVerts
        vi = mesh.v(qi,:);
        ni = mesh.n(qi,:);
        [ot,ov] = onering(mesh, qi);
        nNbrs = numel(ov);
        
        dists = vmag(vadd( mesh.v(ov,:), -vi));
        if ( smooth_radius == -1 )
            r = max(dists) * 1.5;
        else
            r = smooth_radius;
        end

        for k = 1:nNbrs
            t = dists(k) / r;
            p = max(2*t*t*t - t*t + 1, 0);
            ni = ni + p * mesh.n(ov(k),:);
        end
        N(qi,:) = normalize(ni);
        
    end
end



