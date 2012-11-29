% read a test mesh
mesh = readMesh('patch.obj');
nVtx = 63;

mesh = readMesh('patch2.obj', 'n');
nVtx = 414;

plotMesh(mesh, 'vefb');

% number of per-vertex neighbours in dijkstra graph
Kgraph = 8;

% make local neighbour distance matrix (does double-duty as connectivity graph)
X = mesh.v';
[D,N] = size(X);
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
G = spalloc(N, N, Kgraph * N);
[sorted,index] = sort(distance);
for jj = 1:N
    nbrs = index(2:(1+Kgraph),jj);
    nbrdists = sorted(2:(1+Kgraph),jj);
    G(jj,nbrs) = nbrdists;
end


% full dijkstra's algorithm
D = dijkstra(G, nVtx);
D(D~=Inf)
max(D(D~=Inf))

% dijkstra's algorithm that terminates after max_dist is hit
%  ! still processes full input/output arrays - is this more or less
%    efficient than sending smaller arrays to dijkstra() ?? -RMS
max_dist = 0.015;
D = dijkstra_maxdist(G, nVtx, max_dist);
D(D~=Inf)
max(D(D~=Inf))


% same as above, but function also returns which upwind neighbour
% each point was determined from (ie shortest-paths information)
max_dist = Inf;
[D,Nearest] = dijkstra_maxdist(G, nVtx, max_dist);
D(D~=Inf)
max(D(D~=Inf))