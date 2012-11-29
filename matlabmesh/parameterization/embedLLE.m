function [ uv, xyz ] = embedLLE( points, W, iboundarypts, boundaryloop, mesh )
%[ uv, xyz ] = embedLLE( points, W, iboundarypts, boundaryloop, mesh )
%   Detailed explanation goes here

% based on http://www.cs.toronto.edu/~roweis/lle/code/lle.m

% The original LLE algorithm assumed that weights were normalized to
% sum to 1. It *seems* like we can skip normalization by constructing
% the M matrix as L'*L, where L is the laplacian matrix (D-W), and then
% solving the generalized eigenvalue problem Lv = \lambda Dv. 
% This doesn't give the exact same result, though - I believe that in this
% modified formulation the vertices are 'weighted' by D.
do_normalize = 0;

[N,D] = size(points);
d = 2;

% normalize weights (maybe)
if do_normalize
    for v = 1:N
        [i,j] = find(W(v,:));
        W(v, j) = W(v, j) / sum(W(v, j));  
    end
    D = speye(N,N);
else
    D = sum(W,2);
    D = sparse(1:N,1:N,D,N,N); 
end


if exist('boundaryloop','var')
    
    % need to construct 2Nx2N system that includes both X and Y 
    % uv-values. We will use form [X,0; 0,Y]
    Ld = [(D-W),zeros(N,N); zeros(N,N), (D-W)];
    
    if exist('mesh','var')
        A = sparse([],[],[],2*N,2*N,4*numel(iboundarypts));
        for ii = 1:numel(iboundarypts)
            i = iboundarypts(ii);    ix = i;    iy = ix+N;

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
    else
      %   faster version that doesn't need mesh
        A = sparse([],[],[],2*N,2*N,4*numel(boundaryloop));
        loop = boundaryloop;
        for ii = 1:numel(loop)
            jx = loop(ii);
            jy = jx + N;
            kx = loop( mod(ii,numel(loop)) + 1 );
            ky = kx + N;

            A(jx,ky) = A(jx,ky) + 1;
            A(kx,jy) = A(kx,jy) - 1;
            A(jy,kx) = A(jy,kx) - 1;
            A(ky,jx) = A(ky,jx) + 1;
        end
    end

    %M = Ld'*Ld;
    M = (Ld-A)'*(Ld-A);
    %M = 1*(Ld'*Ld) - A;
   

    fudge = 10e-8;      % numerical fudge from sec 3.4
    %fudge = -0.00088;      % numerical fudge from sec 3.4
    
    if ~isposdef(M + fudge*speye(2*N,2*N))
        fprintf('M is not posdef w/ fudge - really bad!\n');
    end    
   
    
    B = sparse([],[],[],2*N,2*N,2*numel(iboundarypts));
    Bi = [iboundarypts, iboundarypts+N];
    B(sub2ind(size(B),Bi,Bi)) = 1;  
    
    % solve eigensystem
    eigsoptions.disp = 0; 
    eigsoptions.isreal = 1; 
    eigsoptions.issym = issymmetric(M);
    neigs = 5;
%    [V,D] = eigs(B, M + fudge*speye(2*N,2*N), neigs, 'lm', eigsoptions);
    [V,D] = eigs(M + fudge*speye(2*N,2*N), B, neigs, 0, eigsoptions);
%    [V,D] = eigs(M + fudge*speye(2*N,2*N), neigs, 0, eigsoptions);
%    [V,D] = eigs(M + fudge*speye(2*N,2*N), neigs, 0, eigsoptions);

    k = 3;
    uv = [V(1:N, k), V(N+1:2*N, k)];  
    uv = real(uv);
    xyz = [];


else

    % [RMS] improve conditioning of weight matrix (from [Mullen08])
    % [RMS] by setting this to pretty big values we can overcome some of the
    %   problems where the solution is 'rotated' in 3D...
    fudge = 10e-8;
    %fudge = 0.007;     % for doghead
    %fudge = 0.03;     % for doghead
    D = D + fudge*speye(N,N);

    % construct matrix M
    Ld = D-W;
    M = Ld' * Ld;
    issym = issymmetric(M);    

    if ~isposdef(M)
        fprintf('M is not posdef - bad!\n');
    end        
    
    if exist('iboundarypts','var')
        B = sparse(iboundarypts,iboundarypts,ones(numel(iboundarypts),1),   N,N,   numel(iboundarypts));
        D = D .* B;
    end    
    
    % CALCULATION OF EMBEDDING
    options.disp = 0; 
    options.isreal = 1; 
    options.issym = issym; 
    [Y,eigenvals] = eigs(M, D, d+2,0,options);

    uv = [ Y(:,2), Y(:,3) ];
    xyz = [ Y(:,2), Y(:,3), Y(:,1) ];
end

end
