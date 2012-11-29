function [ uv, xyz ] = embedLEM( points, W, iboundarypts )
% [ uv, xyz ] = embedLEM(points, W)  
%   Laplacian EigenMaps [Belkin02]
%    points is set of high-dimensional points
%    W is weight matrix
%    uv is 2D embedding coordinates
%    xyz is 3D embedding cooordinates
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]

% instead of taking 2 lowest eigenvectors as (u,v), construct a 2Nx2N
% matrix and take the lowest non-zero eigenvector as (u;v).
% [NOTE] this doesn't work! 
solve_single_row = 0;

[N,dim] = size(points);

% construct matrices
D = sparse([],[],[],N,N,N);
D(1:N+1:N*N) = sum(W,2);  % set D(i,i) = -sum_j W(i,j)


if solve_single_row
    W = [W,sparse(N,N); sparse(N,N), W];
    D = [D,sparse(N,N); sparse(N,N), D];
    L = D - W;
    fudge = 10e-8;      % numerical fudge [Mullen02]
    L = L + fudge*speye(2*N,2*N);
else
    L = D - W;
    fudge = 10e-8;      % numerical fudge [Mullen02]
    L = L + fudge*speye(N,N);
end

issym = issymmetric(L);
if ~issym
    fprintf('[embedLEM] WARNING: Laplacian L is not symmetric - LEM may fail...\n');
elseif ~isposdef(L)
    fprintf('[embedLEM] WARNING: Laplacian L is not posdef - LEM may fail...\n');
end


% solve only for boundary points
if exist('iboundarypts','var')
    B = sparse(iboundarypts,iboundarypts,ones(numel(iboundarypts),1),   N,N,   numel(iboundarypts));
    D = D.*B;
end

% solve eigensystem
eigsoptions.disp = 0; 
eigsoptions.isreal = 1; 
eigsoptions.issym = issym;
d = 2;      % dimension we want to embed into
[V,eigenvals] = eigs(L, D, d+2, 0, eigsoptions);

V = real(V);

if solve_single_row
    k=4;
    uv = [V(1:N, k), V(N+1:2*N, k)]; 
    xyz = points;
else
    e2 = real(eigenvals(2,2));
    e3 = real(eigenvals(3,3));
    uv = real( [V(:,2)/sqrt(e2), V(:,3)/sqrt(e3)] );
    fprintf('*** [embedLEM] scaling eigenvectors...\n');
%    uv = real( [V(:,2), V(:,3)] );
    xyz = [V(:,2), V(:,3), V(:,1)];   % embed in 3 dimensinos
end


end
