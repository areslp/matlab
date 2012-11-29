function [ theta ] = vangle( v1, v2 )
%ANGLE vangle( A, B ) 
%  compute angle between vectors A and B

% [TODO] there are more accurate ways to compute an angle than acos

v1 = normalize(v1);
v2 = normalize(v2);
d = vdot(v1,v2);
d = clamp(d,-1,1);
theta = acos(d);

end
