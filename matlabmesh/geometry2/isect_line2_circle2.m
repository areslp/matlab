function [ hits ] = isect_line2_circle2( p0, p1, c, r )
%ISECT_LINE2_CIRCLE2 compute intersection points between line and circle
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   based on http://www.geometrictools.com/Documentation/IntersectionLine2Circle2.pd


d = normalize(p1 - p0);
k = p0 - c;

discrim = vdot(d,k)^2 - vmag2(d)*(vmag2(k)-r^2);
if discrim < 0
    hits = [];
elseif discrim >= 0
    t0 = (-vdot(d,k) + sqrt(discrim)) / vmag2(d);
    t1 = (-vdot(d,k) - sqrt(discrim)) / vmag2(d);
    hits = [ p0 + t0*d ; p0 + t1*d ];
end
