function [ A ] = vertexArea( mesh, vertices, mode )
%[ A ] = vertexArea( mesh, vertices, mode )
%  estimate area of mesh vertices
%   Valid modes are:
%     'uniform' - area = 1   
%     'onering' - sum of areas of one-ring triangles
%     'voronoi' - voronoi areas from [Meyer02]  (invalid for obtuse triangles)
%     'mixed'   - mixed voronoi/geometric areas from [Meyer02]
%  If vertices empty or undefined, compute areas for entire mesh 

if ~ exist('mode', 'var')
    mode = 'onering';
end

if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = 1:size(mesh.v,1);
end
n = numel(vertices);

if strcmp(mode,'voronoi') || strcmp(mode,'mixed')
    W = cotanWeights(mesh, vertices);
end

% [TODO] move loops inside strcmps (much faster)
if strcmp(mode,'onering')
    A = zeros(n,1);
    for i = 1:n
        A(i) = vertexArea_OneRing( mesh, vertices(i) );
    end
elseif strcmp(mode,'voronoi')
    A = zeros(n,1);
    for i = 1:n
        A(i) = vertexArea_Voronoi( mesh, vertices(i) );
    end    
elseif strcmp(mode,'mixed')    
    A = zeros(n,1);
    for i = 1:n
        A(i) = vertexArea_Mixed( mesh, vertices(i) );
    end      
else
    A = ones(n,1);
end
    

end


% W = vector of cotan weights
function [ A ] = vertexArea_Voronoi( mesh, v, W )
    nbrs = find(W ~= 0);
    A = 0;
    for ni = 1:numel(nbrs)
        j = nbrs(ni);
        A = A + W(j) * vmag2(mesh.v(v,:) - mesh.v(j,:));
    end
    A = A / 8;
end



% W = vector of cotan weights
function [ A ] = vertexArea_Mixed( mesh, v, W )
    A = 0;

    [ot,ov] = onering(mesh,v);
    for ti = 1:numel(ot);
        f = mesh.f(ot(ti),:);
        if triobtuse(mesh.v(f,:))
            vother = f(f~=v);
            vtri = vadd(mesh.v(vother,:), -mesh.v(v,:));
            angle = vangle(vtri(1,:), vtri(2,:));
            if angle > pi/2
                A = A + triarea( mesh.v(f(1),:), mesh.v(f(2),:), mesh.v(f(3),:) ) / 2;
            else
                A = A + triarea( mesh.v(f(1),:), mesh.v(f(2),:), mesh.v(f(3),:) ) / 4;
            end
        else
            % [TODO] replace w/ vcot
            vother = f(f~=v);
            P = mesh.v(v,:);
            Q = mesh.v(vother(1),:);
            R = mesh.v(vother(2),:);
            angleQ = vangle( normalize(P-Q), normalize(R-Q) );
            angleR = vangle( normalize(P-R), normalize(Q-R) );
            A = A + (1/8)*( vmag2(P-R)*cot(angleQ) + vmag2(P-Q)*cot(angleR) );
        end
    end
end



function [ A ] = vertexArea_OneRing( mesh, v )
    A = 0;
    [ot,ov] = onering(mesh,v);
    for i = 1:numel(ot)
        f = mesh.f(ot(i),:);
        A = A + triarea( mesh.v(f(1),:), mesh.v(f(2),:), mesh.v(f(3),:) );
    end
end