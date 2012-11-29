function [ distance ] = implicitPlane( x, y, z, origin, normal )
%IMPLICITPLANE returns signed distance to plane

vec = vadd([x,y,z],-origin);
distance = dot(vec,normal);

end
