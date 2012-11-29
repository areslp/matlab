function [ A, isboundary ] = pointArea( P, vertices, mode, options )
% [ A, isboundary ] = pointArea( P, vertices, mode, options )
% estimate area of points on surface
%   Valid modes are:
%     'uvVoronoi' - voronoi cell area in UV-space
%     'uvDelArea' - full delaunay area in UV-space (doesn't set isboundary)
%     'uvDelAreaR' - delaunay area in radius (doesn't set isboundary)
%                         - needs options.radius parameter
%     'uvDelRing' - delaunay one-ring area in UV-space (ditto)
%        -( above need sparse matrices options.u, options.v)
%    
%  If vertices empty or undefined, compute areas for entire mesh 

if ~ exist('mode', 'var')
    mode = 'uv';
end

if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = 1:size(P,1);
end
n = numel(vertices);

A = zeros(n,1);
isboundary = zeros(n,1);
if strcmp(mode, 'uvVoronoi')
    for i = 1:n
        vi = vertices(i);        
        idx = find(options.u(vi,:));
        x = [0, full(options.u(vi,idx))];
        y = [0, full(options.v(vi,idx))];
        uv = [x',y'];

        %voronoi(x,y);   % plot 2D voronoi diagram
        [v,c] = voronoin([x(:) y(:)]);
        v1 = v(c{1},:);
        
        % voronoi cell area 
        A(vi) = polyarea(v1(:,1),v1(:,2));
        
        % voronoi area is undefined if cell is
        % not closed. In that case, use (one-ring area)/3
        if isnan(A(vi)) 
            A(vi) = 0;
            TRI = delaunay(x,y);
            [rr,cc]=find(TRI==1);
            tris = TRI(rr,:);
            for ti = 1:numel(rr)
                triv = uv(tris(ti,:),:);
                A(vi) = A(vi) + triarea(triv(1,:), triv(2,:), triv(3,:));
            end
            A(vi) = A(vi) / 3;
        end

        % Cells for points on boundary would ideally be connected 
        % to 'infinity' vertex. But sometimes the cell is actually
        % closed, at a point very far away. Using a radius threshold
        % is an estimate for when that happens...
        dists = vmag2(uv);
        r = 1.5*mean(dists);
        dists = vmag2(v1);
        if max(dists) > r
%            voronoi(x,y);
%            pause;
            isboundary(vi) = 1;
        end    
    end
elseif strcmp(mode, 'uvVoronoi3')
    error('not done yet...\n');
    for i = 1:n
        vi = vertices(i);        
        idx = find(options.u(vi,:));
        x = [0, full(options.u(vi,idx))];
        y = [0, full(options.v(vi,idx))];
        uv = [x',y'];

        %voronoi(x,y);   % plot 2D voronoi diagram
        [v,c] = voronoin([x(:) y(:)]);
        v1 = v(c{1},:);
        
        % voronoi cell area 
        A(vi) = polyarea(v1(:,1),v1(:,2));
        
        % voronoi area is undefined if cell is
        % not closed. In that case, use (one-ring area)/3
        if isnan(A(vi)) 
            A(vi) = 0;
            TRI = delaunay(x,y);
            [rr,cc]=find(TRI==1);
            tris = TRI(rr,:);
            for ti = 1:numel(rr)
                triv = P(tris(ti,:),:);
                A(vi) = A(vi) + triarea(triv(1,:), triv(2,:), triv(3,:));
            end
            A(vi) = A(vi) / 3;
        end

        % Cells for points on boundary would ideally be connected 
        % to 'infinity' vertex. But sometimes the cell is actually
        % closed, at a point very far away. Using a radius threshold
        % is an estimate for when that happens...
        dists = vmag2(uv);
        r = 1.5*mean(dists);
        dists = vmag2(v1);
        if max(dists) > r
%            voronoi(x,y);
%            pause;
            isboundary(vi) = 1;
        end    
    end    
elseif strcmp(mode, 'uvDelArea')
    for i = 1:n
        vi = vertices(i);        
        idx = find(options.u(vi,:));
        x = [0, full(options.u(vi,idx))];
        y = [0, full(options.v(vi,idx))];
        uv = [x',y'];
        TRI = delaunay(x,y);
        n = size(TRI,1);
        for ti = 1:n
            triv = uv(TRI(ti,:),:);
            A(vi) = A(vi) + triarea(triv(1,:), triv(2,:), triv(3,:));
        end
    end
elseif strcmp(mode, 'uvDelAreaR')
    for i = 1:n
        vi = vertices(i);        
        idx = find(options.u(vi,:));
        x = [0, full(options.u(vi,idx))];
        y = [0, full(options.v(vi,idx))];
        dist2 = x.^2 + y.^2;
        nidx = 1:numel(x);
        nidx = nidx(dist2 < options.radius*options.radius);
        uv = [x(nidx)',y(nidx)'];
        
        [K, A(vi)] = convhulln(uv); 
        
%         TRI = delaunay(x,y);
%         n = size(TRI,1);
%         for ti = 1:n
%             triv = uv(TRI(ti,:),:);
%             if ( vmag(triv(1,:)) < options.radius/2 && vmag(triv(2,:)) < options.radius/2 && vmag(triv(3,:)) < options.radius/2 )
%                 A(vi) = A(vi) + triarea(triv(1,:), triv(2,:), triv(3,:));
%             end
%         end
    end
elseif strcmp(mode, 'uvDelRing')
   for i = 1:n
        vi = vertices(i);        
        idx = find(options.u(vi,:));
        x = [0, full(options.u(vi,idx))];
        y = [0, full(options.v(vi,idx))];
        N = numel(x);
        uv = [x',y'];
        TRI = delaunay(x,y);
        mesh = makeMesh([uv,zeros(N,1)],TRI, zeros(N,1));
        A(vi) = vertexArea(mesh, 1, 'onering');
    end   
elseif strcmp(mode, 'uvDelRing3')
   for i = 1:n
        vi = vertices(i);        
        idx = find(options.u(vi,:));
        x = [0, full(options.u(vi,idx))];
        y = [0, full(options.v(vi,idx))];
        N = numel(x);
        uv = [x',y'];
        TRI = delaunay(x,y);
        xyz3 = [P(vi,:); P(idx,:)];
        mesh = makeMesh(xyz3,TRI, zeros(N,1));
        A(vi) = vertexArea(mesh, 1, 'onering');
    end   
end


end