function [ distance ] = implicitSphere( x, y, z, center, radius )
%IMPLICITSPHERE returns signed distance to sphere (negative == inside)

dist = vmag(vadd(center, -[x,y,z]));
distance = vadd(dist,-radius);

end
