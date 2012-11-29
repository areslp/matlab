function [ uv, EC, cons1, cons2 ] = embedDNCP( mesh, W, constrain_vtx1, constrain_vtx2, alpha )
% [uv,EC] = embedDNCP(mesh, W, constrain_vtx)  
%   Discrete Natural Conformal Parameterization [Desbrun02]
%     W: optional custom weight matrix (used instead of cotan weights)
%     constrain_vtx1/2: pin these vertices
%     EC: returned conformal energy u'(D-A)u
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]

% list of boundary vertex indices
iboundary = mesh.vidx(mesh.isboundaryv~=0);

if~exist('alpha','var')
    alpha = 1;
end

if ~exist('constrain_vtx1','var') | constrain_vtx1 < 1
    cons1 = iboundary(1);
else
    cons1 = constrain_vtx1;
end

% pick two vertices to constrain. Use vertex 1, and furthest
% vertex from vertex 1 ( [Desbrun02] says to use two furthest
% vertices, but that is O(N^2)! - this is a reasonable approx.)
if ~exist('constrain_vtx2','var') | constrain_vtx2 < 1
    dists = vmag2(vadd(mesh.v(iboundary,:), -mesh.v(cons1,:)));
    [maxval,maxi] = max(dists);
    cons2 = iboundary(maxi);    
else
    cons2 = constrain_vtx2;
end
ifixed = [cons1,cons2];
fixedUV = [cons1, 0,0;  cons2,alpha,0];


% construct initial NxN weight matrix
N = numel(mesh.vidx);
if ~exist('W','var') | isempty(W)
    W = cotanWeights(mesh);
end
W = -W;
W(1:N+1:N*N) = -sum(W,2);  % set W(i,i) = -sum_j W(i,j)




% need to construct 2Nx2N system that includes both X and Y 
% uv-values. We will use form [X,0; 0,Y]
Ld = [W,sparse(N,N); sparse(N,N), W];
rhs = zeros(2*N,1);

% construct natural boundary conditions matrix [Desbrun02]
% [TODO] all but the entries that are between pairs of edge verts cancel
%   out to 0....could maybe do this more efficiently...
% A = sparse([],[],[],2*N,2*N,4*numel(iboundary));
% for ii = 1:numel(iboundary)
%     i = iboundary(ii);    ix = i;    iy = ix+N;
%     
%     vtris = oneringf(mesh,i);
%     for ti = 1:numel(vtris)
%         f = mesh.f(vtris(ti), :);
%         [j,k] = tripick(f,i);
%         jx = j;  jy = j+N;
%         kx = k;  ky = k+N;
%         
%         A(ix,ky) = A(ix,ky) - 1;
%         A(ix,jy) = A(ix,jy) + 1;
% 
%         A(iy,jx) = A(iy,jx) - 1;
%         A(iy,kx) = A(iy,kx) + 1;
%     end
%     
% end


% faster version based on loops (but need to support multiple boundary loops!
A = sparse([],[],[],2*N,2*N,4*numel(iboundary));
for li = 1:numel(mesh.loops)
    loop = mesh.loops{li};
    for ii = 1:numel(loop)
        jx = loop(ii);
        jy = jx + N;
        kx = loop( mod(ii,numel(loop)) + 1 );
        ky = kx + N;

        A(jx,ky) = A(jx,ky) + 1;
    %    A(ky,jx) = A(ky,jx) + 1;
        A(kx,jy) = A(kx,jy) - 1;
    %    A(jy,kx) = A(jy,kx) - 1;    
    end
end
A = (A + A');


% temp
% ux = mesh.v(:,1);
% uy = mesh.v(:,2);
% uv = [ux,uy];
% u = [ux;uy];
% EDEdges = 0;
% [i,j]=find(mesh.e);
% for ii = 1:numel(i);
%     v1 = i(ii);  v2 = j(ii);
%     ti = mesh.te(v1,v2);
%     if ti ~= 0
%         f = mesh.f(ti,:);
%         v3 = f(f~=v1 & f~=v2);
%         e1 = normalize(vadd(mesh.v(v1,:),-mesh.v(v3,:)));
%         e2 = normalize(vadd(mesh.v(v2,:),-mesh.v(v3,:)));
%         len = vmag2( uv(v1,:) - uv(v2,:) );
%         EDEdges = EDEdges + (1/4)*vcot(e1,e2)*len;
%     end
% end
% ED1 = ( ux'*W*ux + uy'*W*uy );
% ED2 = u'*Ld*u;
% A = u'*A*u;
% EC1 = ED1-A;
% EC2 = ED2-A;
% temp


% construct natural conformal system L_C   [Mullen08]
Lc = Ld - A;
LcCons = Lc;

% set fixed-vertex rows to M(i,i) = 1
ifixed = [ifixed,ifixed+N];
LcCons(ifixed,:) = 0;
LcCons(sub2ind(size(Lc),ifixed, ifixed)) = 1;

% set fixed vertex positions
rhs(fixedUV(:,1)) = fixedUV(:,2);
rhs(fixedUV(:,1)+N) = fixedUV(:,3);

% move w_i*u_i column values of constrained verts over to RHS
% (this is not strictly necessary, but it keeps the matrix symmetric...)
rhsadd = zeros(size(rhs));
for k = 1:numel(ifixed)
    ci = ifixed(k);
    col = LcCons(:,ci);
    col(ci) = 0;
    rhsadd = rhsadd + rhs(ci)*col;
end
LcCons(:,ifixed) = 0;
LcCons(sub2ind(size(Lc),ifixed, ifixed)) = 1;
rhs = rhs-rhsadd;

% solve system and unpack uv-values
uv = LcCons \ rhs;

EC = 0.5*uv'*Lc*uv;
%fprintf(' E_D: %f   A: %f  E_C: %f  E_C+cons: %f\n', uv'*Ld*uv,  uv'*A*uv, EC, uv'*LcCons*uv);

uv = [uv(1:N), uv(N+1:2*N)];

% return this as energy...not continuous, but smooth between big jumps...
EC = log(EC.^2) + turningNumber(uv(mesh.loops{1},:));

end
