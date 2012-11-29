function [ ringv, vpos, ringt ] = oneringExpMap( mesh, nVert )
%ONERINGEXPMAP Summary of this function goes here
%   Detailed explanation goes here

[t,v] = onering(mesh, nVert, 'ccw');

vcenter = mesh.v(nVert,:);
e1 = vadd(mesh.v(v,:),-vcenter);
e2 = vadd(mesh.v(circshift(v,[0,-1]),:),-vcenter);
angles = vangle(e1,e2);
scale = (2*pi) / sum(angles);
vangles = cumsum(angles * scale );
vangles = vangles-vangles(1);
vr = vmag(e1);
vpos = [vr.*cos(vangles),vr.*sin(vangles),zeros(numel(v),1)];
ringv = v;
ringt = t;

end
