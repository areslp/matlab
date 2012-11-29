function [ deformed_mesh ] = fair_Laplacian( mesh )
%FAIR_LAPLACIAN solve for zero-length laplacian vectors

deformed_mesh = mesh;

Bi = mesh.vidx( mesh.isboundaryv == 1 );  % boundary verts
wbdry = 100;
constraints = [Bi, mesh.v(Bi,:), repmat(wbdry, numel(Bi), 1) ];


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
%D = [ M*mesh.v(:,1), M*mesh.v(:,2), M*mesh.v(:,3) ];
D = zeros(size(mesh.v));

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
deformed_mesh.n = estimateNormal(deformed_mesh);




