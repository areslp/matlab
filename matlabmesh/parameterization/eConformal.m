function [ EConformal, EDirichlet, Area, Turning ] = eConformal( surface, W )
% [ EConformal, EDirichlet, Area, Turning ] = eConformal( surface, W )

N = numel(surface.vidx);
W = -W;
W(1:N+1:N*N) = -sum(W,2);  % set W(i,i) = -sum_j W(i,j)
Ld = [W,sparse(N,N); sparse(N,N), W];

% set boundary conditions
A = sparse(2*N,2*N);
nloops = numel(surface.loops,1);
for k = 1:nloops
    loop = surface.loops{k};
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

uv = [surface.u(:,1); surface.u(:,2)];

EDirichlet = 0.25*uv'*Ld*uv;
Area = 0.25*uv'*A*uv;
EConformal = EDirichlet - Area;
Turning = turningNumber(surface.u(surface.loops{1},:));

end
