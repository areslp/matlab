function [ deformed_mesh, rotations ] = deformRotInvCoords( mesh, consF, consP )
%[ deformed_mesh ] = deformRotInvCoords( mesh, consF, consP )
%
%   consF:      rows of [vtx_i, m11, m12,... m33, weight_i]   (m = 3x3 transformation matrix)
%   consP:      rows of [vtx_i, x, y, z, weight_i]

pcount = numel(mesh.vidx);

[A,B,N] = meshFrames(mesh);             % [A,B,N] = tangent vectors and normals (Nx3 lists)
[C1,C2,C3] = edgeCoeffs(mesh, A,B,N);   % [alpha,beta,rho] from survey eq.(23)  (stored as adjacency-structured sparse matrices)

% set up linear system for rotations

% construct linear system M.V = [rhsA,rhsB,rhsN]
% M is (edges*3)x(points*3) system matrix
% start col for vert i = (i-1)*3+1;
[ii,jj] = find(mesh.e);
ecount = numel(ii);
Mrot = sparse([],[],[],  3*ecount, 3*pcount, 3*ecount*4);
RHSrot = zeros(3*ecount, 3);
matidx = 1;
for k = 1:ecount
    i = ii(k);  j = jj(k);
    ci = (i-1)*3+1;  cj = (j-1)*3+1;
    rij = matidx;
    matidx = matidx+3;
 
    % construct row (Fi-Fj) - A*Fi = (1-A)Fi - Fj = Rij*Fi - Fj = 0   (are the signs right?)
    
    Fi = [ A(i,:); B(i,:); N(i,:) ]';
    Fj = [ A(j,:); B(j,:); N(j,:) ]';
    Rij = Fj'*Fi;    

    Mrot(rij  ,ci:ci+2) = Rij(1,:);
    Mrot(rij+1,ci:ci+2) = Rij(2,:);
    Mrot(rij+2,ci:ci+2) = Rij(3,:);

    Mrot(rij,  cj) = -1; 
    Mrot(rij+1,cj+1) = -1; 
    Mrot(rij+2,cj+2) = -1;
end


% compute residual to test system
% testX = zeros(3*pcount,3);
% for k = 1:pcount
%     pi = (k-1)*3+1;
%     testX(pi:pi+2,1:3) = [A(k,:); B(k,:); N(k,:)];
% end
% fprintf('residual is %f\n', sum(sum(abs( Mrot*testX - RHSrot))) );

% construct normal equations
RHSrot = Mrot'*RHSrot;
Mrot = Mrot'*Mrot;

% append soft constraints to system
nCons = size(consF,1);
consi = zeros(3*nCons,1);
consv = zeros(3*nCons,3);
consw = zeros(3*nCons,1);
for i = 1:nCons
    xi = consF(i,1);
    wi = consF(i,11);
    rotation = reshape( consF(i,2:10), 3, 3 )';
    F = [A(xi,:); B(xi,:); N(xi,:)];
    Fcons = mvmul(rotation,F);
    crow = (i-1)*3+1;
    mrow = (xi-1)*3+1;
    consi(crow:crow+2) = [mrow:mrow+2];
    consv(crow:crow+2,:) = Fcons;
    consw(crow:crow+2) = wi;
end
[Mcons,RHScons] = softConstrain(Mrot,RHSrot,consi,consv,consw);

% solve for frames
Frames = Mcons \ RHScons;

% unpack into A2,B2,N2 and ortho-normalize frames
A2 = zeros(size(A));  B2 = A2;  N2 = A2; 
for xi = 1:pcount
    ci = (xi-1)*3+1;
    a = Frames(ci,1:3);    b = Frames(ci+1,1:3);    n = Frames(ci+2,1:3);
    a = normalize(a); b = normalize(b); n = normalize(n);

    % orthonormalize frames (preserving n here...should probably distribute error equally
    a = ncross(b,n);
    b = ncross(n,a);
    
    A2(xi,:) = a;    B2(xi,:) = b;    N2(xi,:) = n;
end



rotations = zeros(pcount,9);
for i = 1:pcount
    Fi = [ A(i,:); B(i,:); N(i,:) ];
    Fj = [ A2(i,:); B2(i,:); N2(i,:) ];
    Rij = (Fj'*Fi);
    rotations(i,:) = reshape(Rij,1,9);
end


% set up linear system for positions
% [TODO] this would be way more efficient if we construct [ii,jj,ss] vectors and passed them directly to sparse()
[ii,jj] = find(mesh.e);
ecount = numel(ii);
Mpos = sparse([],[],[], ecount,pcount,2*ecount);
RHSpos = zeros(ecount,3);
for k = 1:ecount
    i = ii(k);  j = jj(k);
    Mpos(k,i) = -1;
    Mpos(k,j) = 1;
    rhsv = C1(i,j)*A2(i,:) + C2(i,j)*B2(i,:) + C3(i,j)*N2(i,:);
    RHSpos(k,:) = rhsv;
end

% construct normal equations
% Mpos ends up being 2*Laplacian, so we can use a cotan laplacian instead
%  (not sure why we need the factor of two, though...)
RHSpos = Mpos'*RHSpos;
Mpos = Mpos'*Mpos;
%Mpos = -makeOneRingWeights(mesh,'cotan',0)*2;     % why *2 ??
%Mpos(1:pcount+1:pcount*pcount) = -sum(Mpos,2); 

% append soft constraints to system
nCons = size(consP,1);
consi = consP(:,1);
consv = consP(:,2:4);
consw = consP(:,5);
[Mcons,RHScons] = softConstrain(Mpos, RHSpos, consi,consv,consw);

% solve for positions
Xpos = Mcons \ RHScons;

deformed_mesh = mesh;
deformed_mesh.v = Xpos(1:pcount,:);

% set normals to estimated frame normals (useful for debugging)
deformed_mesh.n = N2;

end





% compute mesh frame vectors (tan1,tan2,normal) = (A,B,N) and orientation bits O
function [ A, B, N ] = meshFrames( mesh )

A = zeros(size(mesh.v)); B = A; N = A;
O = zeros(size(mesh.vidx));

%vAreas = faceArea(mesh);
%vNormals = faceNormal(mesh);

nv = numel(mesh.vidx);
for i = 1:nv

    % use mesh normal
    normal = mesh.n(i,:);
    
    % compute area-weighted average normal
    %vTris = oneringf(mesh,i);
    %normals = vNormals(vTris,:) .* repmat(vAreas(vTris,:),1,3);
    %normal = normalize( sum(normals) );
    
    nbrs = find(mesh.e(i,:)>0);
    k = nbrs(1);
    xik = mesh.v(k,:) - mesh.v(i,:);
    xikbar = cross( normal, cross(xik, normal) );

    N(i,:) = normal;
    A(i,:) = normalize(xikbar);
    B(i,:) = cross( normal, A(i,:) );
end

end



% [alpha,beta,rho] from survey eq.(23)  (stored as adjacency-structured sparse matrices)
function [alpha,beta,rho] = edgeCoeffs(mesh, A,B,N)

alpha = mesh.e; beta = alpha; rho = alpha;

nv = numel(mesh.vidx);
for i = 1:nv
    nbrj = find(mesh.e(i,:)>0);
    nbrv = vadd( mesh.v(nbrj,:), -mesh.v(i,:) );
    
    alpha(i,nbrj) = vdot(nbrv, A(i,:));
    beta(i,nbrj) = vdot(nbrv, B(i,:));
    rho(i,nbrj) = vdot(nbrv, N(i,:));
end

end


