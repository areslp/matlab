function [ x, y, z ] = sph2cartY( THETA, PHI, R )
% [ x, y, z ] = sph2cartY( THETA, PHI, R )
%   consistent w/ mathworld, graphics papers, etc
x = R .* cos(THETA) .* sin(PHI);
y = R .* sin(THETA) .* sin(PHI);
z = R .* cos(PHI);
end

