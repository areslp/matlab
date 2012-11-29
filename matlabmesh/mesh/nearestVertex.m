function [ index, distance ] = nearestVertex( mesh, position )
%NEARESTVERTEX find nearest vertex to position

dists = vmag2( vadd(mesh.v, -position) );
[distance,index] = min(dists);

end
