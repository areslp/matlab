function [ Tpq ] = sphereNormCoords( radius, p, q, tan1, tan2 )
% [ u, v ] = sphereNormCoords( radius, p, q, tan1, tan2 )
%   returns normal coordinates (u,v) of q in tangent space at p on sphere
%   of given radius. Basis of TS is [tan1, tan2, normalize(p)]  (in 3D)

if vmag(p-q) < sqrt(eps)
    Tpq = [0,0];
    return;
end

% find great-circle
plane_norm = ncross( normalize(p), normalize(q) );
[e1,e2] = tangentFrame(plane_norm);
pc = [ vdot(e1,p), vdot(e2,p) ];
qc = [ vdot(e1,q), vdot(e2,q) ];
r_g = vangle(pc,qc) * radius;

v = normalize(p-q);
Tv = [vdot(v,tan1), vdot(v,tan2)];
theta_g = atan2( Tv(2), Tv(1) );

Tpq = [ r_g * cos(theta_g), r_g * sin(theta_g) ];


end
