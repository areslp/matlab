function [ deformed_mesh ] = meshDeform_Laplacian( mesh, constraints, target_normals )
%MESHDEFORM_LAPLACIAN basic Laplacian mesh deformation 
%
%   constraints:      rows of [vtx_i, x, y, z, weight_i]
%   target_normals:   pre-determined normal field over desired target surface

deformed_mesh = mesh;


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


% solve for deformed mesh (need local nbrhoods to get tangents below)
deformed_mesh.v(:,1) = M \ D(:,1);
deformed_mesh.v(:,2) = M \ D(:,2);
deformed_mesh.v(:,3) = M \ D(:,3);


% transform D vectors
if exist('target_normals','var')
    
    orig_D = D;    
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

        % construct new frame at deformed vertex, with target normal
        pi = deformed_mesh.v(qi,:);
        ni = target_normals(qi,:);
        pj = deformed_mesh.v(qj,:);
        uij = vadd(pj,-pi);
        uij = vadd(uij, -vdot(uij,ni)*ni);
        uij = normalize(uij);
        e2 = ncross(ni, uij);   

        new_d = alpha*ni + beta*uij + gamma*e2;

        D(qi,:) = new_d;
    end
    
    % solve again with new D vectors
    deformed_mesh.v(:,1) = M \ D(:,1);
    deformed_mesh.v(:,2) = M \ D(:,2);
    deformed_mesh.v(:,3) = M \ D(:,3);
    
end

deformed_mesh.n = estimateNormal(deformed_mesh);




