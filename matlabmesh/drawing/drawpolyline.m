function [ output_args ] = drawpolyline( polyline, linewidth, mode )
%DRAWPOLYLINE draw a polyline
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
% mode specifies flags for rendering (eg 've')
%   'e' - draw edges
%   'v' - draw vertices

if ~ exist('linewidth')
    linewidth = 1;
end

if ~exist('mode','var')
    mode = 'e';
end

[n,dim] = size(polyline);

if findstr(mode,'v')
    if dim == 2
        scatter( polyline(:,1), polyline(:,2), 40, 'filled' );    
    elseif dim == 3
        scatter3( polyline(:,1), polyline(:,2), polyline(:,3), 40, 'filled' );    
    end
end

if findstr(mode, 'e')
    if dim == 2
        line(polyline(:,1),polyline(:,2), 'LineWidth', linewidth);
    elseif dim == 3
        line(polyline(:,1),polyline(:,2), polyline(:,3), 'LineWidth', linewidth);
    end
end

end
