function [ vi ] = pickVertex( mesh, redraw )
%PICKVERTEX interactively pick a vertex in mesh
%   returns vertex index

if ~ exist('redraw','var') | redraw ~= 0
    plotMesh(mesh, 'vefb');
end

fig = gcf;
datacursormode on;
waitforbuttonpress;
dcm_obj = datacursormode(fig);
f = getCursorInfo(dcm_obj);
datacursormode off;
vpos = f.Position;
vi = nearestVertex(mesh, vpos);

end
