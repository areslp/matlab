function [ centers ] = circle2_2pr( p1, p2, r )
%CIRCLE2_2PR find circle centers from 2 points and radius (2 results)
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   based on http://mathforum.org/library/drmath/view/53027.html

q = sqrt( sum( (p2-p1).^2 ) );
h = (p1+p2)/2;
m = [ p1(2)-p2(2), p2(1)-p1(1) ];
d = sqrt(r^2-(q/2)^2) * m / q;
centers = [h+d;h-d];

end
