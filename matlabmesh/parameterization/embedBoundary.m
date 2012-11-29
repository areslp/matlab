function [ boundaryUV ] = embedBoundary( mesh, mode )
%  [ boundaryUV ] = embedBoundary( mesh, mode )
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   - valid modes are:
%     'circle' - arc-length embedding on circle
%     'copy'   - current UV values

boundary = mesh.loops{1};
boundaryUV = [];

% find perimeter length
n = size(boundary,1);
k1 = [1:n];
k2 = [2:n,1];
vdiff = mesh.v(boundary(k1),:) - mesh.v(boundary(k2),:);
lengths = sqrt( sum( vdiff .* vdiff, 2) );
dists = cumsum(lengths);

% fixed translation and scale
offset = [0,0];
scale = 1;

if strcmp(mode, 'circle')
    r = dists(n) / (2*pi);
    r = r * scale;
    udists = dists / dists(n);
    thetas = udists * 2 * pi;
    boundaryUV = [boundary, r*cos(thetas) + offset(1), r*sin(thetas) + offset(2)];
elseif strcmp(mode, 'square')

    %error('square mode is not implemented yet\n');

    % first map to circle
    r = dists(n) / (2*pi) * scale;
    thetas = (dists / dists(n)) * 2 * pi;
    
    edge_origins = [ 0,0; 1,0; 1,1; 0,1 ];
    edge_dirs = [ 1,0; 0,1; -1,0; 0,-1 ];
    
    for k = 0:3
        idx = find(thetas > k*pi/2 & thetas <= (k+1)*pi/2);
        d = dists(idx);

        if k < 3
            next_idx = find(thetas > (k+1)*pi/2 & thetas <= (k+2)*pi/2);
            next_d = dists(next_idx);
        else
            next_d = dists(1) + dists(end);
        end
        
        d = (d-d(1)) / (next_d(1)-d(1));
        d
        origin = edge_origins(k+1,:);
        dir = edge_dirs(k+1,:);
        edgeUV = [boundary(idx), origin(1) + scale*d*dir(1) + offset(1), origin(2) + scale*d*dir(2) + offset(2) ];
        boundaryUV = [boundaryUV ; edgeUV]; 
    end        

    
elseif strcmp(mode, 'copy')
    boundaryUV = [boundary, mesh.u(boundary,:)];
end


