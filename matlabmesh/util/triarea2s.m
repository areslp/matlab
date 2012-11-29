function [ A ] = triarea2s( p1, p2, p3 )
% [ A ] = triarea2s( p1, p2, p3 )
%   signed area of triangle defined by 3 2D points
    A = 0.5 * ( (p2(1)-p1(1))*(p3(2)-p1(2)) - (p2(2)-p1(2))*(p3(1)-p1(1)) );
end
