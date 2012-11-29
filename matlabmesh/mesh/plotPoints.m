function [ ] = plotPoints( pointset, mode, vC )
% [ ] = plotPoints( pointset, mode, vC )
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2009]
% mode specifies flags for rendering (eg 'v')
%   'v' = draw vertex points
%   'n' = draw normals
%   'u' = use uv-values as vertices
%   'U' = use uv-values as vertices, and try to consistently-orient/scale plots
%   'i' = plot vertex indices (numbers)  ** WARNING: very slow 
%   'O' = pretty-plot for output (thicker lines, etc)
%   vC is vertex or face colors (passed as patch() FaceVertexCData option)
%      - can be set to rgb triplet or scalar value, in which case the
%        current colormap is used to map to rgb

% hardcoded pointset color (should be a parameter...)
faceColorRGB = [1,0,0];
edgeColor = 'black';
boundaryColor = 'blue';
lineWidth = 1;
borderWidth = 3;

if ~exist('mode','var')
    mode = 'v';
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

v = pointset.v;
n = pointset.n;
nVerts = size(pointset.v,1);

if findstr(mode, 'u')
   v = [ pointset.u, zeros(nVerts,1) ];
   minx = min(v(:,1));  maxx = max(v(:,1));
   miny = min(v(:,2));  maxy = max(v(:,2));
   pad = max(maxx-minx,maxy-miny)/15;
   axis([minx-pad maxx+pad miny-pad maxy+pad]);
   mode = [mode,'v'];
elseif findstr(mode,'U')           % try to keep parameterization plot consistently oriented...
   align = normalize(pointset.u(1,:) - pointset.u(2,:));
   theta = atan2(align(2),align(1));
   rotate = [cos(theta),sin(theta); -sin(theta),cos(theta)];
   v = [ mvmul(rotate, pointset.u), zeros(nVerts,1) ];
   center = mean(v);
   v = vadd(v,-center);
   minx = min(v(:,1));  maxx = max(v(:,1));
   miny = min(v(:,2));  maxy = max(v(:,2));
   vwidth = max(v(:,1)) - min(v(:,1));
   vheight = max(v(:,2)) - min(v(:,2));
   vscale = max(vwidth,vheight);
%   v = v/vscale;
%   axis([-1 1 -1 1]);
   axis([1.1*minx 1.1*maxx 1.1*miny 1.1*maxy]);
   mode = [mode,'v'];
end


if findstr(mode, 'v') 
    lightdir = normalize([0.5,0.5,0.5]);
    if ~isempty(n)
        NdotL = abs(vdot(n, lightdir));
        ptColors = [NdotL,NdotL,NdotL];
    else
        ptColors = [1,0,0]
    end
    scatter3(v(:,1),  v(:,2), v(:,3), 100, ptColors, 'Marker','.');
end


if findstr(mode, 'l')
    pos = pointset.bounds(3,:) + 5 * (pointset.bounds(2,:)-pointset.bounds(3,:));
    %light('Position',pos,'Style','infinite');
    for az = [0,90,180,270]
        h = light('Style','infinite');
        lightangle(h, az,45)
    end
    lighting gouraud;
end


%if findstr(mode, 'e')
%    patch('Vertices', v, 'Faces', pointset.f, 'FaceColor', 'none', 'EdgeColor', edgeColor, 'EdgeLighting', 'none', 'LineWidth', lineWidth);
%end


%if findstr(mode, 'b')
%   nLoops = size(pointset.loops,2); 
%   for i = 1:nLoops
%       vi = [pointset.loops(:,i); pointset.loops(1,i)];
%       line( v(vi,1), v(vi,2), v(vi,3), 'LineWidth', borderWidth, 'Color', boundaryColor );
%   end
%end


if findstr(mode, 'n') & size(pointset.n,1) == nVerts
    nlen = 0.05 * vmag( pointset.bounds(1,:) - pointset.bounds(2,:) );
    for i = 1:nVerts
        vi = [pointset.v(i,:); pointset.v(i,:) + (nlen * pointset.n(i,:))];
        line( vi(:,1), vi(:,2), vi(:,3), 'LineWidth', 1, 'Color', 'green' );
    end
end       

if findstr(mode,'i')
    nlen = 0.025 * vmag( pointset.bounds(1,:) - pointset.bounds(2,:) );
    for i = 1:nVerts
        vi = [v(i,:) + (nlen * pointset.n(i,:)); v(i,:)];
        text( vi(1,1), vi(1,2), vi(1,3), int2str(pointset.vidx(i)) );
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
