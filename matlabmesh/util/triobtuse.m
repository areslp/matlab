function [ isobtuse ] = triobtuse( tripts )
%TRIOBTUSE return 1 if any triangle angles are obtuse
%   tripts = [3xdim] vector of triangle points

% [TODO]  http://mathworld.wolfram.com/ObtuseTriangle.html describes a
%    possibly more efficient obtuse test, based on edge-lengths...

i = [1,2,3];
ip = [3,1,2];
in = [2,3,1];

vp = normalize(tripts(ip,:) - tripts(i,:));
vn = normalize(tripts(in,:) - tripts(i,:));
angles = vangle(vp,vn);

isobtuse = ~isempty(angles(angles>pi/2));

end
