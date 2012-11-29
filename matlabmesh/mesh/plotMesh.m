function [ ] = plotMesh( mesh, mode, vC )
% [ ] = plotMesh( mesh, mode, vC )
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
% mode specifies flags for rendering (eg 'veb')
%   'v' = draw vertex points
%   'e' = draw triangle edges
%   'b' = draw boundary edges
%   'n' = draw normals
%   'f' = draw filled faces
%   'l' = apply lighting to faces
%   'u' = use uv-values as vertices
%   'U' = use uv-values as vertices, and try to consistently-orient/scale plots
%   'i' = plot vertex indices (numbers)  ** WARNING: very slow 
%   'O' = pretty-plot for output (thicker lines, etc)
%   vC is vertex or face colors (passed as patch() FaceVertexCData option)
%      - can be set to rgb triplet or scalar value, in which case the
%        current colormap is used to map to rgb

% hardcoded mesh color (should be a parameter...)
faceColorRGB = [1,0,0];
edgeColor = 'black';
boundaryColor = 'blue';
lineWidth = 1;
borderWidth = 3;

if ~exist('mode','var')
    mode = 'efb';
end

newplot;
bWasHeld = ishold;
hold all;

if findstr(mode,'O')
    lineWidth = 1;
    borderWidth = 3;
    % nice vals for screenshots when boundary is important to see
    %lineWidth = 2;
    %borderWidth = 15;
    faceColorRGB = [0.7,0.7,0.7];
end

v = mesh.v;
n = mesh.n;
nVerts = size(mesh.v,1);

if findstr(mode, 'u')
   n = zeros(nVerts,3);
   v = [ mesh.u, zeros(nVerts,1) ];
   minx = min(v(:,1));  maxx = max(v(:,1));
   miny = min(v(:,2));  maxy = max(v(:,2));
   pad = max(maxx-minx,maxy-miny)/15;
   axis([minx-pad maxx+pad miny-pad maxy+pad]);
elseif findstr(mode,'U')           % try to keep parameterization plot consistently oriented...
   n = zeros(nVerts,3);
   align = normalize(mesh.u(1,:) - mesh.u(2,:));
   theta = atan2(align(2),align(1));
   rotate = [cos(theta),sin(theta); -sin(theta),cos(theta)];
   v = [ mvmul(rotate, mesh.u), zeros(nVerts,1) ];
   %center = mean(v)
   center = [ 0.5*(min(v) + max(v)) ];
   v = vadd(v,-center);
   minx = min(v(:,1));  maxx = max(v(:,1));
   miny = min(v(:,2));  maxy = max(v(:,2));
   vwidth = max(v(:,1)) - min(v(:,1));
   vheight = max(v(:,2)) - min(v(:,2));
   vscale = max(vwidth,vheight);
%   v = v/vscale;
%   axis([-1 1 -1 1]);
   axis([1.1*minx 1.1*maxx 1.1*miny 1.1*maxy]);
end
   
if findstr(mode, 'v')
    offset = 0.01;
    voff = v + offset * n;
    scatter3( voff(:,1), voff(:,2), voff(:,3), 40, 'filled' );
end


if findstr(mode, 'l')
    pos = mesh.bounds(3,:) + 5 * (mesh.bounds(2,:)-mesh.bounds(3,:));
    %light('Position',pos,'Style','infinite');
    for az = [0,90,180,270]
        h = light('Style','infinite');
        lightangle(h, az,45)
    end
    lighting gouraud;
end

if findstr(mode, 'f')
    if exist('vC')
        if size(vC,2) == 3   % truecolor colors
            if size(vC,1) == size(v,1)
                patch('Vertices', v, 'VertexNormals', n, 'Faces', mesh.f, 'FaceVertexCData', vC, 'FaceColor','interp', 'EdgeColor', 'none');
            else
                patch('Vertices', v, 'VertexNormals', n, 'Faces', mesh.f, 'FaceVertexCData', vC, 'FaceColor', 'flat', 'EdgeColor', 'none');
            end
        else
            colorbar;
            if size(vC,1) == size(v,1)
                patch('Vertices', v, 'VertexNormals', n, 'Faces', mesh.f, 'FaceVertexCData', vC, 'FaceColor','interp', 'EdgeColor', 'none');
            else
                patch('Vertices', v, 'VertexNormals', n, 'Faces', mesh.f, 'FaceVertexCData', vC, 'FaceColor', 'flat', 'EdgeColor', 'none');
            end
        end
    else
        vC = repmat(faceColorRGB, size(mesh.v,1), 1);
        patch('Vertices', v, 'VertexNormals', n, 'Faces', mesh.f, 'FaceVertexCData', vC, 'FaceColor','interp', 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'BackFaceLighting', 'lit', 'AmbientStrength',0.01, 'DiffuseStrength',0.4, 'SpecularStrength',0.4);
    end
end


if findstr(mode, 'e')
    patch('Vertices', v, 'Faces', mesh.f, 'FaceColor', 'none', 'EdgeColor', edgeColor, 'EdgeLighting', 'none', 'LineWidth', lineWidth);
end


if findstr(mode, 'b')
   nLoops = numel(mesh.loops); 
   for i = 1:nLoops
       vi = mesh.loops{i}; 
       vi = [vi(1:end); vi(1)];
       line( v(vi,1), v(vi,2), v(vi,3), 'LineWidth', borderWidth, 'Color', boundaryColor );
   end
end



if findstr(mode, 'n') & size(mesh.n,1) == nVerts
    nlen = 0.05 * vmag( mesh.bounds(1,:) - mesh.bounds(2,:) );
    for i = 1:nVerts
        vi = [mesh.v(i,:); mesh.v(i,:) + (nlen * mesh.n(i,:))];
        line( vi(:,1), vi(:,2), vi(:,3), 'LineWidth', 1, 'Color', 'green' );
    end
end       

if findstr(mode,'i')
    nlen = 0.025 * vmag( mesh.bounds(1,:) - mesh.bounds(2,:) );
    for i = 1:nVerts
%        if ~mesh.isboundaryv(i)
%            continue;
%        end
        vi = [v(i,:) + (nlen * mesh.n(i,:)); v(i,:)];
        text( vi(1,1), vi(1,2), vi(1,3), int2str(mesh.vidx(i)) );
        line( vi(:,1), vi(:,2), vi(:,3), 'LineWidth', 1, 'Color', 'green' );
    end   
end
   
if ~findstr(mode, 'u') & ~findstr(mode, 'U')
    view([45 35])
    axis vis3d;
else
    axis equal;
end

if ~bWasHeld
    hold off;   
end
