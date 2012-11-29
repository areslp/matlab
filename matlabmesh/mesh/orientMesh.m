function [ outmesh ] = orientMesh( mesh, mode )
%ORIENTCCW Summary of this function goes here
%   'ccw2' - mesh is in xy plane, use z=1 as normal

outmesh = mesh;
n = size(outmesh.f,1);

if strcmp(mode, 'ccw2')
    for ti = 1:n
       t = outmesh.f(ti,:);
       v = outmesh.v(t,:);
       area = cross2( v(2,:)-v(1,:), v(3,:)-v(1,:) );
       if area < 0
           outmesh.f(ti,:) = [t(2), t(1), t(3)];
       end
    end
    
end