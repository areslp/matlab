function [ THETA, PHI, R ] = cart2sphY( x, y, z )
% [ THETA, PHI, R ] = cart2sphY( x, y, z )
%   consistent w/ mathworld, graphics papers, etc
R = sqrt( x.^2 + y.^2 + z.^2 );
THETA = atan2( y, x );
PHI = acos( z ./ R );
end
