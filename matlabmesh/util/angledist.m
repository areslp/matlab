function [ anglediff ] = angledist( a1, a2 )
%ANGLEDIST Summary of this function goes here
%   Detailed explanation goes here

v1 = [cos(a1),sin(a1)];
v2 = [cos(a2),sin(a2)];
anglediff = vangle(v1,v2);

end
