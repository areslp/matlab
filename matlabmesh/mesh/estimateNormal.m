function [ N ] = estimateNormal( surface, vertices, mode )
%[ N ] = estimateNormal( mesh, vertices, mode )
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   - if vertices is undefined/empty, compute for entire mesh
%   - valides modes are:
%       'faceavg' - straight average of one-ring face normals
%       'faceavg_area' - area-weighted one-ring face normals

if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = 1:size(surface.v,1);
end
if ~ exist('mode','var')
    mode = [];
end

if isfield(surface,'f') & size(surface.f,1) > 0
    N = estimateNormal_mesh(surface,vertices,mode);
else
    N = estimateNormal_pointset(surface,vertices,mode);
end

end







%% estimate normals for point set
function [N] = estimateNormal_pointset( pointset, vertices, mode)


n = numel(vertices);
dim = size(pointset.v,2);
   
N = zeros(n,dim);
for i = 1:n
    nbrs = find(pointset.e(i,:));
    Qip = [pointset.v(i,:);pointset.v(nbrs,:)];
    qmean = mean(Qip);
    [U,E,V]=svd(vadd(Qip,-qmean));
end

end








%% estimate normals for mesh
function [N] = estimateNormal_mesh( mesh, vertices, mode)

if ~ exist('mode', 'var') | isempty(mode)
    mode = 'faceavg';
end


n = numel(vertices);
dim = size(mesh.v,2);
   
N = zeros(n,dim);
if strcmp(mode, 'faceavg')
    for i = 1:n
        N(i,:) = estimateNormal_faceAvg( mesh, vertices(i) );
    end
elseif strcmp(mode, 'faceavg_area')
    for i = 1:n
        N(i,:) = estimateNormal_faceAreaAvg( mesh, vertices(i) );
    end    
elseif strcmp(mode, 'nbrsvd')
    for i = 1:n
        N(i,:) = estimateNormal_SVD( mesh, vertices(i) );
    end
end

end



function [ N ] = estimateNormal_faceAvg( mesh, v )
    [ot,ov] = onering(mesh, v);
    sumn = [0,0,0];
    for k = 1:numel(ot)
        f = mesh.f(ot(k),:);
        e1 = mesh.v(f(2),:) - mesh.v(f(1),:);
        e2 = mesh.v(f(3),:) - mesh.v(f(1),:);
        sumn = sumn + ncross( e1, e2 );
    end
    N = normalize(sumn);
end


function [ N ] = estimateNormal_faceAreaAvg( mesh, v )
    [ot,ov] = onering(mesh, v);
    sumn = [0,0,0];
    for k = 1:numel(ot)
        f = mesh.f(ot(k),:);
        e1 = mesh.v(f(2),:) - mesh.v(f(1),:);
        e2 = mesh.v(f(3),:) - mesh.v(f(1),:);
        kross = cross(e1,e2);
        sumn = sumn + normalize(kross) * 0.5 * vmag( kross );
    end
    N = normalize(sumn);
end




function [ N ] = estimateNormal_SVD( mesh, i )
    nbrs = find(mesh.e(i,:));
    Qip = [mesh.v(i,:);mesh.v(nbrs,:)];
    qmean = mean(Qip);
    [U,E,V]=svd(vadd(Qip,-qmean));
    N = V(3,:);
    N = normalize(N);
end
