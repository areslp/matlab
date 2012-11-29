function [ pointset ] = makePointSet( v, n, u, c )
%[ pointset ] = makePointSet( v, n, u )
% initialize pointset data structure from vertex/normal/uv-arrays
%
% pointset.v = Nx3 list of vertices
% pointset.n = Nx3 list of vertex normals
% pointset.u = Nx2 list of vertex uv's
% pointset.vidx = Nx1 list of vertex indices (useful for various things)
% pointset.bounds = [min;max;center] bounding box
% pointset.e = symmetric sparse array of edges
%             pointset.e(i,j) = |v_j-v_i| if connected, 0 otherwise
% pointset.valence = valence of each vertex of pointset
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2009]

pointset = struct('v', [], 'vidx', [], 'n', [], 'u', [], 'c', [], 'e', [], 'bounds', []);

% pointset version string (useful for checking types of cached pointsetes)
pointset.version = globalConfig('pointsetVersion');

pointset.v = v;
nVerts = size(pointset.v,1);

if exist('n','var') && size(n,1) == nVerts
	pointset.n = n;
end
if exist('u','var') && size(u,1) == nVerts
	pointset.u = u;
end
if exist('c','var') && size(c,1) == nVerts
	pointset.c = c;
end

pointset.vidx = (1:nVerts)';


% compute pointset bounds -  min, max, center
pointset.bounds = [ min(pointset.v) ; max(pointset.v) ; 0.5*(min(pointset.v) + max(pointset.v)) ];


% construct K-nn edge matrix
K = 8;
E = sparse([],[],[],nVerts,nVerts,K*nVerts);
for i = 1:nVerts
    distances = vmag(vadd(pointset.v,-pointset.v(i,:)));
    [sorted,index] = sort(distances);
    nbrs = index(2:(K+1));
    E(i,nbrs) = sorted(2:(K+1));
end
pointset.e = E;
    
% compute valences (and one-rings, etc?)
pointset.valence = sum(E>0,2);

end
