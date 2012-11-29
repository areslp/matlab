function [ THETA, PHI, R ] = cart2sphM( x, y, z )
% [ THETA, PHI, R ] = cart2sphM( x, y, z )
%   consistent w/ mathematica Spherical coordinate system
%     x(theta,phi) = cos(phi)*sin(theta)
%     y(theta,phi) = sin(phi)*sin(theta)
%     z(theta,phi) = cos(theta)
R = sqrt( x.^2 + y.^2 + z.^2 );
PHI = atan2( y, x );
THETA = acos( z ./ R );
end
