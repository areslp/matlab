function [ uv ] = embedSCP( mesh, mode, W, faceAreaWeighted )
% uv = embedSCP(mesh, mode, W)  
%   Spectral Conformal Parameterization [Mullen08]
%    mode: computation mode
%      'fiedler'     - Fielder vector          (sec 3.2)
%      'generalized' - generalized (default)   (sec 3.3)
%      'robust'      - generalized w/ search for isposdef(Lc+eps*I) matrix
%      'lle'         - LLE embedding
%      'test'        - test embedding
%    W: weight matrix (will compute cotan weights if not specified)
%    faceAreaWeighted: if =1, use face-area-weighted Area matrix
%    
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]

if sum(mesh.isboundaryv) == 0
    uv = embedSSCP(mesh);
    return;
end


if ~exist('mode','var')
    mode = 'generalized';
end
if ~exist('faceAreaWeighted', 'var')
    faceAreaWeighted = 0;
end

if ~exist('W','var') | isempty(W)
    W = cotanWeights(mesh);
end
W = -W;   % negate so that sums on diagonal are positive

% list of boundary vertex indices
iboundary = mesh.vidx(mesh.isboundaryv~=0);

% construct initial NxN weight matrix
N = numel(mesh.vidx);
W(1:N+1:N*N) = -sum(W,2);  % set W(i,i) = -sum_j W(i,j)

% need to construct 2Nx2N system that includes both X and Y 
% uv-values. We will use form [X,0; 0,Y]
Ld = [W,sparse(N,N); sparse(N,N), W];

% construct natural boundary conditions matrix [Desbrun02]
% [TODO] all but the entries that are between pairs of edge verts cancel
%   out to 0....could maybe do this more efficiently...
if ~faceAreaWeighted
    A = sparse([],[],[],2*N,2*N,4*numel(iboundary));
    for ii = 1:numel(iboundary)
        i = iboundary(ii);    ix = i;    iy = ix+N;

        vtris = oneringf(mesh,i);
        for ti = 1:numel(vtris)
            f = mesh.f(vtris(ti), :);
            [j,k] = tripick(f,i);
            jx = j;  jy = j+N;
            kx = k;  ky = k+N;

            A(ix,ky) = A(ix,ky) - 1;
            A(ix,jy) = A(ix,jy) + 1;

            A(iy,jx) = A(iy,jx) - 1;
            A(iy,kx) = A(iy,kx) + 1;
        end
    end
    

%     % faster version based on loops (but need to support multiple boundary loops!
%     A2 = sparse([],[],[],2*N,2*N,4*numel(iboundary));
%     loop = mesh.loops;
%     for ii = 1:numel(loop)
%         jx = loop(ii);
%         jy = jx + N;
%         kx = loop( mod(ii,numel(loop)) + 1 );
%         ky = kx + N;
% 
%         A2(jx,ky) = A2(jx,ky) + 1;
%         A2(ky,jx) = A2(ky,jx) + 1;
%         A2(kx,jy) = A2(kx,jy) - 1;
%         A2(jy,kx) = A2(jy,kx) - 1;
%     end    
    
else
    A = Ld;      A(:) = 0;
    areas = faceArea(mesh);
    
    [i,j] = find(mesh.e);
    for k = 1:numel(i)
        sum1 = 0;
        sum2 = 0;
        f1 = mesh.te(i(k),j(k));
        if f1 ~= 0
            Tijk = areas(f1);
            sum1 = sum1 + 1/Tijk;
            sum2 = sum2 + -1/Tijk;
        end
        f2 = mesh.te(j(k),i(k));
        if f2 ~= 0
            Tijl = areas(f2);
            sum1 = sum1 + -1/Tijl;
            sum2 = sum2 + 1/Tijl;
        end            
        ui = i(k);  vi = ui+N;
        uj = j(k);  vj = uj+N;
        A(ui,vj) = (1/2)*sum1;
        A(vi,uj) = (1/2)*sum2;
    end
    A = 2*A;   % need to balance Ld....
end



% construct natural conformal system L_C   [Mullen08]
Lc = Ld-A;

if strcmp(mode, 'fiedler')
    % solve eigensystem
    fudge = 10e-8;      % numerical fudge from sec 3.4
    eigsoptions.disp = 0; 
    eigsoptions.isreal = 1; 
    eigsoptions.issym = issymmetric(Lc);
    neigs = 5;
    [V,D] = eigs(Lc + fudge*speye(2*N,2*N), neigs, 0, eigsoptions);
    k = 3;
    uv = [V(1:N, k), V(N+1:2*N, k)];
    
elseif strcmp(mode, 'generalized')

    Bi = [mesh.loops{1}, mesh.loops{1}+N];
    B = sparse(Bi,Bi,ones(numel(Bi),1),2*N,2*N);
%    B = sparse([],[],[],2*N,2*N,2*numel(iboundary));
%    Bi = [iboundary, iboundary+N];
%    B(sub2ind(size(B),Bi,Bi)) = 1;

    fudge = 10e-8;      % numerical fudge from sec 3.4
%    if ~isposdef(Lc + fudge*speye(2*N,2*N))
       % fprintf('Lc is not posdef w/ fudge - bad!\n');
%    end    
    
    % solve eigensystem
    eigsoptions.disp = 0; 
    eigsoptions.isreal = 1; 
    eigsoptions.issym = issymmetric(Lc);
    neigs = 5;
    [V,D] = eigs(Lc + fudge*speye(2*N,2*N), B, neigs, 0, eigsoptions);
%    [V,D] = eigs(B, Lc + fudge*speye(2*N,2*N), neigs, 'lm', eigsoptions);
    k = 3;
    uv = real([V(1:N, k), V(N+1:2*N, k)]);  
%    x1 = real(V(1:N, k));   y1 = imag(V(1:N, k));
%    x2 = real(V(N+1:2*N, k));   y2 = imag(V(N+1:2*N, k));
%    xo = real(V(1:N, k));  yo = real(V(N+1:2*N, k));
%    uv = [x1,y1];
    
elseif strcmp(mode, 'robust')
    B = sparse([],[],[],2*N,2*N,2*numel(iboundary));
    Bi = [iboundary, iboundary+N];
    B(sub2ind(size(B),Bi,Bi)) = 1;

    eigsoptions.disp = 0; 
    eigsoptions.isreal = 1; 
    eigsoptions.issym = issymmetric(B);    
    
    % solve eigensystem
    k = 7;
    fudge = 10^-k;
    LcR = Lc+fudge*speye(2*N,2*N);
    while ~isposdef(LcR) & k > 1
        k = k-1;
        fudge = 10^-k;
        LcR = Lc+fudge*speye(2*N,2*N);
    end

    neigs = 5;
    [V,D] = eigs(B, LcR, neigs, 'lm', eigsoptions);
    idx = find(D>0);
    k = idx(1);
%k=4;
    uv = [V(1:N, k), V(N+1:2*N, k)];    
    
elseif strcmp(mode, 'lle')
    B = sparse([],[],[],2*N,2*N,2*numel(iboundary));
    Bi = [iboundary, iboundary+N];
    B(sub2ind(size(B),Bi,Bi)) = 1;
    
    fudge = 10e-3;      
    LcR = Lc + fudge*speye(2*N,2*N);  % numerical fudge from sec 3.4  
    M = LcR' * LcR;
    issym = issymmetric(M);
    
    options.disp = 0; 
    options.isreal = 1; 
    options.issym = issym; 
%    [V,D] = eigs(M, 5, 0, options);
    [V,D] = eigs(M, B, 5, 0, options);

    k = 3;
    uv = [V(1:N, k), V(N+1:2*N, k)];     
    
elseif strcmp(mode,'test')
    uv = embedSSCP(mesh);
end

end

















% 3D embedding 
%   This does 3D Laplacian Eigenmaps right now. Since we assume
%   mesh is closed (no boundary), x/y/z components are disconnected
%   anyway. But what if we had a boundary? Then we could include 
%   mesh area A matrix (would not be same as above though, since it
%   is 3D). What would happen then? And what if we added constraints?
function [ uv ] = embedSSCP( mesh )

W = -cotanWeights(mesh);

% construct initial NxN weight matrix
N = numel(mesh.vidx);
W(1:N+1:N*N) = -sum(W,2);  % set W(i,i) = -sum_j W(i,j)

% solve eigensystem
fudge = 10e-8;      % numerical fudge from sec 3.4
eigsoptions.disp = 0; 
eigsoptions.isreal = 1; 
eigsoptions.issym = issymmetric(W);
neigs = 5;
[V,D] = eigs(W + fudge*speye(N,N), neigs, 0, eigsoptions);

k = 2;
uv = [V(1:N,k), V(1:N,k+1), V(1:N,k+2)];

end
