function [ uv ] = embedPSCP( pointset, W )
% [ uv ] = embedPSCP( pointset, W )
%   Pointset extension of Spectral Conformal Parameterization [Mullen08]
%    W: weight matrix 
%    
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]


W = -W;   % negate so that sums on diagonal are positive

% construct initial NxN weight matrix
N = numel(pointset.vidx);
W(1:N+1:N*N) = -sum(W,2);  % set W(i,i) = -sum_j W(i,j)

% need to construct 2Nx2N system that includes both X and Y 
% uv-values. We will use form [X,0; 0,Y]
Ld = [W,sparse(N,N); sparse(N,N), W];

% construct natural boundary conditions matrix [Desbrun02]
for li = 1:numel(pointset.loops)
    loop = pointset.loops{1};
    A = sparse([],[],[],2*N,2*N,4*numel(loop));
    for ii = 1:numel(loop)
        jx = loop(ii);
        jy = jx + N;
        kx = loop( mod(ii,numel(loop)) + 1 );
        ky = kx + N;

        A(jx,ky) = A(jx,ky) + 1;
        A(ky,jx) = A(ky,jx) + 1;
        A(kx,jy) = A(kx,jy) - 1;
        A(jy,kx) = A(jy,kx) - 1;
    end    
end 

% construct natural conformal system L_C   [Mullen08]
Lc = Ld-A;

    

Bi = [loop, loop+N];
%B = sparse([],[],[],2*N,2*N,2*numel(loop));
B = sparse(Bi,Bi,ones(numel(Bi),1),2*N,2*N);
%B(sub2ind(size(B),Bi,Bi)) = 1;

fudge = 10e-8;      % numerical fudge from sec 3.4
%if ~isposdef(Lc + fudge*speye(2*N,2*N))
   % fprintf('Lc is not posdef w/ fudge - bad!\n');
%end    

% solve eigensystem
eigsoptions.disp = 0; 
eigsoptions.isreal = 1; 
eigsoptions.issym = issymmetric(Lc);
neigs = 5;
[V,D] = eigs(Lc + fudge*speye(2*N,2*N), B, neigs, 0, eigsoptions);
%[V,D] = eigs(B, Lc + fudge*speye(2*N,2*N), neigs, 'lm', eigsoptions);
k = 3;
uv = real([V(1:N, k), V(N+1:2*N, k)]);  

end












