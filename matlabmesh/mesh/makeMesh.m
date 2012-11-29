function [ mesh ] = makeMesh( v, f, n, u, options )
%[ mesh ] = makeMesh( v, f, n, u )
% initialize mesh data structure from vertex/face/normal/uv-arrays
%
% mesh.v = Nx3 list of vertices
% mesh.n = Nx3 list of vertex normals
% mesh.u = Nx2 list of vertex uv's
% mesh.f = Mx3 list of triangles
% mesh.fn = Mx3 list of triangle face normals
% mesh.vidx = Nx1 list of vertex indices (useful for various things)
% mesh.fidx = Nx1 list of face indices (useful for various things)
% mesh.bounds = [min;max;center] bounding box
% mesh.e = symmetric sparse array of edges
%             mesh.e(i,j) = 1 if boundary edge, 2 if interior edge
% mesh.te = sparse array of edge triangles
%             mesh.e(i,j) = tri1, mesh.e(j,i) = tri2
%             for boundary tri, one of these will be 0...
% mesh.loops = boundary loops of mesh
% mesh.valence = valence of each vertex of mesh
% mesh.isboundaryv = 0/1 flags for boundary vertices
% mesh.isboundaryt = 0/1 flags for boundary triangles
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]

mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);

if ~ exist('options', 'var')
    options = struct('noboundary',0); 
end
if ~ isfield(options,'noboundary')
     options.noboundary = 0; 
end

% mesh version string (useful for checking types of cached meshes)
mesh.version = globalConfig('meshVersion');

mesh.v = v;
nVerts = size(mesh.v,1);

if exist('n','var') && size(n,1) == nVerts
	mesh.n = n;
end
if exist('u','var') && size(u,1) == nVerts
	mesh.u = u;
end

mesh.f = f;
nFaces = size(mesh.f,1);


mesh.vidx = (1:nVerts)';
mesh.fidx = (1:nFaces)';


% compute mesh bounds -  min, max, center
mesh.bounds = [ min(mesh.v) ; max(mesh.v) ; 0.5*(min(mesh.v) + max(mesh.v)) ];

% compute face normals
% [RMS TODO] does this give outward-pointing normals??
mesh.fn = zeros(nFaces,3);
for i = 1:nFaces
    fv = mesh.v( mesh.f(i,:), : );
    e1 = fv(2,:) - fv(1,:);
    e2 = fv(3,:) - fv(1,:);
    mesh.fn(i,:) = normalize( cross(e1, e2) );
end

% construct sparse edge matrix
e = zeros(nFaces*3*2,3);
for i = 1:nFaces
    for j = 0:2
        t1 = [ mesh.f(i,j+1), mesh.f(i,mod(j+1,3)+1), 1 ];
        t2 = [ mesh.f(i,mod(j+1,3)+1), mesh.f(i,j+1), 1 ];
%        e = cat(1, e, t1);
%        e = cat(1, e, t2);
        e(((i-1)*6)+(j*2)+1,:) = t1;
        e(((i-1)*6)+(j*2)+2,:) = t2;
    end
end
mesh.e = sparse( e(:,1), e(:,2), e(:,3), nVerts, nVerts );


% ok, elements of mesh.e are 2 if they are interior edges
% and 1 if the are boundary edges. Make sorted, unique list 
% of boundary edge vertex pairs
[i,j] = find(mesh.e==1);
mesh.isboundaryv = zeros(nVerts,1);
mesh.isboundaryv(i) = 1;
be = sortrows( sort([i,j],2) );
be = unique(be, 'rows');


mesh.isboundaryf = ( mesh.isboundaryv(mesh.f(:,1)) + mesh.isboundaryv(mesh.f(:,2)) + mesh.isboundaryv(mesh.f(:,3)) );
mesh.isboundaryf = mesh.isboundaryf > 0;


loops = [];
if options.noboundary == 0
    % ok now construct ordered boundary loops
    % [TODO] use half-edge structures to do this, will
    %  automatically deal with orientation problems
    loopk = 1;
    bloop = [];
    while numel(be) > 0
        bloop = [];
        a1 = be(1,1); a2 = be(1,2);
        be(1,:) = [];
        bloop = cat(1,bloop,a2);
        while size(be,1) > 0
           nextrow = find( be(:,1) == a2 & be(:,2) ~= a1 );
           if nextrow
              b2 = be(nextrow,1); b3 = be(nextrow,2);
           else
              nextrow = find( be(:,2) == a2 & be(:,1) ~= a1 ); 
              b3 = be(nextrow,1); b2 = be(nextrow,2);
           end
           if isempty(nextrow)
               loops{loopk} = bloop;
               loopk = loopk+1;
               break;
           else
               be(nextrow,:) = [];
               bloop = cat(1,bloop,b3);
               a1 = b2; a2 = b3;
           end
        end
    end
    if ~isempty(bloop) 
        loops{loopk} = bloop;
        loopk = loopk+1;
    end
end


% loop construction may have reversed orientation - if so, fix
% [TODO] fix loop construction so this doesn't happen...
for k = 1:numel(loops)
    loop = loops{k};
    prev_idx = [3,1,2];
    loop1 = loop(1);   loop2 = loop(2);
    [fi,fj] = find(mesh.f == loop1);
    for i = 1:numel(fi)
        jp = prev_idx( fj(i) );
        if mesh.f(fi(i), jp) == loop2
            nL = size(loop,1);
            loop = loop(nL:-1:1);
        end
    end   
    loops{k} = loop;
end


% sort loops by size
if ~ isempty(loops)
    for k = 1:numel(loops)
        loopsize(k) = numel(loops{k});
    end
    [sorted,idx] = sort(loopsize,'descend');
    for k = 1:numel(idx)
        mesh.loops{k} = loops{idx(k)};
    end
else
    mesh.loops = [];
end



% now make winged-edge matrix
mesh.te = mesh.e;
mesh.te( find(mesh.e~=0) ) = 0;
for ti = 1:nFaces
    for k = 0:2
        v1 = mesh.f(ti,k+1);  v2 = mesh.f(ti, mod((k+1),3)+1 );
        if mesh.te( v1, v2 ) ~= 0
            tmp = v1; v1 = v2; v2 = tmp;
        end
        mesh.te( v1, v2 ) = ti;
    end
end



% compute valences (and one-rings, etc?)
mesh.valence = zeros(nVerts,1);
for vi = 1:nVerts
   [i,j] = find( mesh.e(vi,:) );    
   mesh.valence(vi) = size(j,2);
end




% do sanity checks
if size(mesh.n,1) == 0
    fprintf('[makeMesh WARNING] no normals\n');
elseif size(mesh.n,1) ~= size(mesh.v,1)
    fprintf('[makeMesh WARNING] # of normals != # of vertices\n');
end
if min(sum(mesh.e)) == 2
    fprintf('[makeMesh WARNING] \"ear\" vertices of valence 2\n');
elseif min(sum(mesh.e)) < 2
    fprintf('[makeMesh WARNING] non-manifold vertices of valence < 2\n');
end
if max(max(mesh.f)) > size(mesh.v,1)
    fprintf('[makeMesh WARNING] faces which reference non-existent vertices\n');
end









%% [RMS] construct half-edge mesh (disabled)
if 0


% mesh.hei is list of linear indices into mesh.he sparse array
%  (convert back to i,j with ind2sub(size(mesh.e), idx)
[i,j] = find(mesh.e ~= 0);
mesh.hei = sub2ind(size(mesh.e),i,j);
heii = 1:numel(mesh.hei);

% mesh.he is sparse array of indices into mesh.hei list
mesh.he = mesh.e;
mesh.he(mesh.hei) = heii;

% hev is 'from' vertex for each halfedge
% henext is next halfedge
mesh.hev = zeros(size(mesh.hei));
mesh.veh = zeros(size(mesh.vidx));
mesh.feh = zeros(size(mesh.fidx));
mesh.henext = zeros(size(mesh.hei));
for fi = 1:size(mesh.f,1);
    f = mesh.f(fi,:);
    f2 = circshift(f,[0,-1]);
    ei = sub2ind(size(mesh.e),f,f2);  % linear indices of face edges
    eii = mesh.he(ei);

    % if any of edges are already assigned, we have to reverse face
    if sum(mesh.hev(eii)) ~= 0
        f = [f(2),f(1),f(3)];
        f2 = circshift(f,[0,-1]);       
        ei = sub2ind(size(mesh.e),f,f2);
        eii = mesh.he(ei);
    end
    
    if sum(mesh.hev(eii)) ~= 0
        error('cannot consistently orient faces...');
    end
    
    mesh.hev(eii) = f;
    mesh.veh(f) = eii;
    mesh.feh(fi) = eii(1);
    mesh.henext(eii) = circshift(eii, [0,-1]);
end

% heflip is flip (reverse) halfedge (0 for boundary edges)
mesh.heflip = mesh.he( sub2ind(size(mesh.e),j,i) );

% ok now have to clean up halfedges that should be 0 because they are boundary edges...
bi = heii(mesh.hev==0);
mesh.hev(bi) = 0;
mesh.henext(bi) = 0;
mesh.heflip(bi) = 0;
mesh.he(mesh.hei(bi)) = 0;
mesh.hei(bi) = 0;

% now set 


end

end
